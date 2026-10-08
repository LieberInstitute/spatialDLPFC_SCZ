#!/usr/bin/env bash
## synthetic integration fixture; never writes to production analysis outputs.
set -eo pipefail
module load conda_R/4.5
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
analysis_dir=${1:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}
Rscript --vanilla - "$analysis_dir" <<'RSCRIPT'
base <- commandArgs(trailingOnly = TRUE)[1]
suppressPackageStartupMessages(library(SpatialExperiment))
source(file.path(base, "config.R"))
source(file.path(base, "R", "core.R"))
source(file.path(base, "R", "stages.R"))
source(file.path(base, "R", "overlap.R"))
main <- function() {
  set.seed(7919)
  root <- tempfile("synthetic_mediation_"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cfg <- config; cfg$min_spots <- 2L; cfg$spots <- "raw.rds"
  env <- list(root = root, out = file.path(root, "outputs"), signature = "synthetic_hitpath_v1")
  d <- expand.grid(SpD = c("spd02", "spd03"), donor = paste0("donor", 1:20), stringsAsFactors = FALSE)
  donor_index <- as.integer(sub("donor", "", d$donor))
  d$Dx <- ifelse(donor_index <= 10, "NTC", "SCZ")
  d$age <- runif(20, 30, 70)[donor_index]
  d$sex <- sample(rep(c("F", "M"), 10))[donor_index]
  d$slide_id <- rep(c("A", "B", "C", "D"), 5)[donor_index]
  d$rin <- runif(20, 6, 9)[donor_index]
  d$section <- paste0("section_", d$donor); d$key <- paste(d$donor, d$SpD, sep = "|")
  spot_group <- rep(seq_len(nrow(d)), each = 5L)
  kind <- rep(1:5, nrow(d))
  cd <- S4Vectors::DataFrame(sample_id = d$section[spot_group], fnl_spd = d$SpD[spot_group],
    vasc_pos = kind %in% c(1, 2, 5), neuropil_pos = kind %in% c(3, 4, 5), neun_pos = kind %in% c(3, 4, 5))
  genes <- data.frame(gene_id = paste0("g", 1:100), gene_name = paste0("Gene", 1:100))
  latent1 <- rnorm(20)[donor_index] + .8 * (d$Dx == "SCZ")
  latent2 <- rnorm(20)[donor_index] + .6 * (d$Dx == "SCZ")
  mu <- matrix(100, 100, length(spot_group))
  mu[1, ] <- 100 * exp(.6 * latent1[spot_group])
  mu[3, ] <- 100 * exp(.6 * latent2[spot_group])
  mu[2, ] <- 100 * exp(.3 * latent1[spot_group] + .2 * latent2[spot_group])
  counts <- matrix(rnbinom(length(mu), mu = mu, size = 40), nrow = 100)
  rownames(counts) <- genes$gene_id; colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
  raw <- SpatialExperiment::SpatialExperiment(assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(genes), colData = cd)
  saveRDS(raw, file.path(root, cfg$spots))
  objects <- lapply(c("vasc", "neuropil", "neun"), function(context) {
    keep <- cd[[paste0(context, "_pos")]]
    z <- Matrix::sparseMatrix(i = seq_len(sum(keep)), j = spot_group[keep], x = 1,
                              dims = c(sum(keep), nrow(d)))
    pb <- as.matrix(counts[, keep] %*% z)
    rownames(pb) <- genes$gene_id; colnames(pb) <- d$key
    ds <- d; ds$ncells <- 3L
    list(counts = pb, logcounts = edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(pb)), log = TRUE),
         samples = ds, genes = genes)
  })
  names(objects) <- c("vasc", "neuropil", "neun")
  h <- data.frame(gene_id = genes$gene_id, p_value_scz = .01, fdr_scz = .1, logFC_scz = .2)
  inputs <- list(objects = objects, historical = setNames(rep(list(h), 3), names(objects)))
  screens <- data.frame(screen_id = c("fgf1_neuropil_vasc", "fgf2_neuropil_vasc"), source = "neuropil",
    target = "vasc", mediator_id = c("g1", "g3"), mediator_symbol = c("Synthetic1", "Synthetic2"), expected_direction = 1)
  tables <- list()
  for (i in 1:2) {
    a <- run_screen(screens[i, ], inputs, cfg, env)
    save_atomic(a, file.path(env$out, "primary", paste0(screens$screen_id[i], ".rds")))
    tables[[i]] <- a$results
  }
  primary <- do.call(rbind, tables)
  ## force only the control-flow trigger in a temporary synthetic fixture.
  ## these labels are not statistical findings and never enter production results.
  primary$screen_hit <- primary$gene_id == "g2"
  write_tsv(primary, file.path(env$out, "primary", "all_pairs.tsv.gz"))
  tryCatch(sensitivity_stage(inputs, screens, cfg, env, workers = 1L), error = function(e) {
    print(read_tsv(file.path(env$out, "sensitivity", "status.tsv")))
    stop(e)
  })
  status <- read_tsv(file.path(env$out, "sensitivity", "status.tsv"))
  assert(all(status$status == "complete"), "Synthetic hit sensitivity failed")
  assert(file.exists(file.path(env$out, "diagnostics", "neuropil_FGF1_FGF2_joint.tsv.gz")), "Joint FGF branch not exercised")
  for (id in screens$screen_id) {
    loo <- read_tsv(file.path(env$out, "diagnostics", paste0(id, "_leave_donor_out.tsv")))
    assert(nrow(loo) == 20 && all(loo$status == "complete"), "Whole-donor influence fixture failed")
    bw <- read_tsv(file.path(env$out, "diagnostics", paste0(id, "_between_within.tsv.gz")))
    assert(nrow(bw) == 100 && sum(bw$primary_hit) == 1, "Between/within full-gene family failed")
  }
  overlap_stage(inputs, screens, cfg, env)
  status <- read_tsv(file.path(env$out, "overlap", "status.tsv"))
  assert(all(status$status == "complete"), "Shared-spot reaggregation fixture failed")
  audit <- read_tsv(file.path(env$out, "overlap", "spot_overlap.tsv"))
  assert(all(audit$shared == 1) && all(audit$n_a == 3) && all(audit$n_b == 3), "Incorrect shared-spot counts")
  for (id in screens$screen_id) for (context in c("vasc", "neuropil")) {
    keep <- read_tsv(file.path(env$out, "overlap", paste0(id, "_", context, "_retention.tsv")))
    assert(all(keep$ncells == 2) && all(keep$retained), "Shared spots not removed from both contexts")
  }
  ## too few remaining spots must be untestable, not a negative association result.
  loss_cfg <- cfg; loss_cfg$min_spots <- 3L
  overlap_stage(inputs, screens, loss_cfg, env)
  loss_status <- read_tsv(file.path(env$out, "overlap", "status.tsv"))
  assert(all(loss_status$status == "untestable_or_failed"), "Insufficient spots misclassified as no hits")
  assert(all(grepl("No pseudobulks remain", loss_status$message)), "Missing sample-loss explanation")
  primary$screen_hit <- FALSE
  write_tsv(primary, file.path(env$out, "primary", "all_pairs.tsv.gz"))
  overlap_stage(inputs, screens, cfg, env)
  status <- read_tsv(file.path(env$out, "overlap", "status.tsv"))
  assert(all(status$status == "not_required_no_primary_hits"), "Zero-hit overlap control flow failed")
  ## an inconsistent raw/pseudobulk input version must fail before reaggregation.
  bad_inputs <- inputs
  bad_inputs$objects$vasc$counts[1, 1] <- bad_inputs$objects$vasc$counts[1, 1] + 1
  mismatch <- try(overlap_stage(bad_inputs, screens, cfg, env), silent = TRUE)
  assert(inherits(mismatch, "try-error") && grepl("does not reproduce", as.character(mismatch)),
         "Raw-count mismatch was not rejected")
  cat("Synthetic positive-trigger diagnostics, joint FGF, raw reconstruction, shared-spot removal, sample-loss, count-mismatch, and zero-hit checks passed.\n")
}
main()
RSCRIPT
