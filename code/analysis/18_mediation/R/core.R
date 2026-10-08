## shared input contracts and inference; no analysis is executed on sourcing.
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
}
read_tsv <- function(path) as.data.frame(data.table::fread(path))
hash <- function(x) digest::digest(x, algo = "sha256")
save_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(x, tmp)
  assert(file.rename(tmp, path), paste("Cannot save", path))
}
canonical_gene <- function(x) sub("\\.[0-9]+$", "", as.character(x))
canonical_spd <- function(x) {
  x <- tolower(as.character(x))
  assert(all(grepl("^spd0[1-7]", x)), "Unrecognized SpD labels")
  substr(x, 1, 5)
}
unique_keys <- function(x, label) {
  assert(!anyNA(x) && all(nzchar(x)) && !anyDuplicated(x), paste(label, "must be unique and nonmissing"))
}
normalize_meta <- function(cd, donor_meta) {
  needed <- c("brnum", "sample_id", "dx", "age", "sex", "ncells")
  assert(all(needed %in% names(cd)), "Missing pseudobulk metadata columns")
  spd_col <- if ("registration_variable" %in% names(cd)) "registration_variable" else "fnl_spd"
  d <- data.frame(donor = as.character(cd$brnum), section = as.character(cd$sample_id),
                  Dx = toupper(as.character(cd$dx)), age = as.numeric(as.character(cd$age)),
                  sex = as.character(cd$sex), SpD = canonical_spd(cd[[spd_col]]),
                  ncells = as.numeric(cd$ncells), stringsAsFactors = FALSE)
  assert(all(d$Dx %in% c("NTC", "SCZ")), "Unexpected diagnosis labels")
  d$key <- paste(d$donor, d$SpD, sep = "|")
  unique_keys(d$key, "Donor-SpD keys")
  assert(all(vapply(split(d$section, d$donor), function(x) length(unique(x)) == 1L, logical(1))),
         "Multiple sections per donor require an explicit aggregation rule")
  ii <- match(d$donor, as.character(donor_meta$brnum))
  assert(!anyNA(ii), "Donor absent from donor_meta")
  for (field in c("slide_id", "rin")) {
    fallback <- donor_meta[[field]][ii]
    if (field %in% names(cd)) {
      v <- cd[[field]]
      both <- !is.na(v) & !is.na(fallback)
      equal <- if (field == "rin") abs(as.numeric(v[both]) - as.numeric(fallback[both])) < 1e-6 else as.character(v[both]) == as.character(fallback[both])
      assert(all(equal), paste("Conflicting", field, "between pseudobulk and donor metadata"))
      v[is.na(v)] <- fallback[is.na(v)]
    } else v <- fallback
    d[[field]] <- if (field == "rin") as.numeric(as.character(v)) else as.character(v)
  }
  for (field in c("dx", "sex", "age")) {
    if (!field %in% names(donor_meta)) next
    v <- donor_meta[[field]][ii]
    observed <- d[[if (field == "dx") "Dx" else field]]
    if (field == "dx") v <- toupper(as.character(v))
    assert(all(is.na(v) | as.character(v) == as.character(observed)), paste("Conflicting donor", field))
  }
  d
}
load_inputs <- function(cfg, root) {
  dm <- read_tsv(file.path(root, cfg$donor_meta))
  unique_keys(dm$brnum, "Donor metadata IDs")
  objects <- lapply(names(cfg$pb), function(context) {
    x <- readRDS(file.path(root, cfg$pb[[context]]))
    rd <- as.data.frame(SummarizedExperiment::rowData(x))
    cd <- as.data.frame(SummarizedExperiment::colData(x))
    genes <- data.frame(gene_id = canonical_gene(rd$gene_id), gene_name = as.character(rd$gene_name))
    unique_keys(genes$gene_id, paste(context, "gene IDs"))
    samples <- normalize_meta(cd, dm)
    counts <- as.matrix(SummarizedExperiment::assay(x, "counts"))
    logcounts <- as.matrix(SummarizedExperiment::assay(x, "logcounts"))
    assert(all(is.finite(counts)) && all(counts >= 0) && all(colSums(counts) > 0), "Invalid counts or empty library")
    rownames(counts) <- rownames(logcounts) <- genes$gene_id
    colnames(counts) <- colnames(logcounts) <- samples$key
    list(counts = counts, logcounts = logcounts, genes = genes, samples = samples)
  })
  names(objects) <- names(cfg$pb)
  historical <- lapply(cfg$historical, function(path) {
    h <- read_tsv(file.path(root, path))
    assert(all(c("ensembl", "p_value_scz", "fdr_scz", "logFC_scz") %in% names(h)), "Historical DEG schema mismatch")
    h$gene_id <- canonical_gene(h$ensembl)
    unique_keys(h$gene_id, "Historical gene IDs")
    h
  })
  list(objects = objects, historical = historical)
}
make_design <- function(d, cfg, extra = character(), include_spd = TRUE) {
  assert(all(complete.cases(d[, c("donor", "Dx", "age", "sex", "SpD", extra), drop = FALSE])), "Missing model metadata")
  donors <- unique(d[, c("donor", "Dx")])
  assert(!anyDuplicated(donors$donor), "Diagnosis differs within donor")
  n <- table(factor(donors$Dx, levels = c("NTC", "SCZ")))
  assert(nrow(donors) >= cfg$min_donors && all(n >= cfg$min_per_dx), "Insufficient independent donors per diagnosis")
  d$Dx <- factor(d$Dx, levels = c("NTC", "SCZ"))
  d$sex <- factor(d$sex)
  d$SpD <- factor(d$SpD)
  if ("slide_id" %in% extra) d$slide_id <- factor(d$slide_id)
  terms <- c("Dx", "age", "sex", if (include_spd) "SpD", extra)
  des <- model.matrix(reformulate(terms, intercept = FALSE), d)
  assert(qr(des)$rank == ncol(des), "Rank-deficient model")
  assert(nrow(des) - ncol(des) >= cfg$min_residual_df, "Insufficient residual degrees of freedom")
  rownames(des) <- d$key
  des
}
match_samples <- function(source, target, cfg, spds = cfg$spds) {
  a <- source$samples; b <- target$samples
  common <- sort(intersect(a$key[a$SpD %in% spds], b$key[b$SpD %in% spds]))
  assert(length(common) > 0, "No matched donor-SpD observations")
  ia <- match(common, a$key); ib <- match(common, b$key)
  for (field in c("donor", "section", "Dx", "age", "sex", "SpD", "slide_id", "rin")) {
    assert(identical(as.character(a[[field]][ia]), as.character(b[[field]][ib])), paste("Source/target disagree on", field))
  }
  d <- b[ib, , drop = FALSE]
  d$source_ncells <- a$ncells[ia]
  d$target_ncells <- b$ncells[ib]
  valid <- complete.cases(d[, c("Dx", "age", "sex", "slide_id", "rin")]) &
    d$source_ncells >= cfg$min_spots & d$target_ncells >= cfg$min_spots
  exclusions <- rbind(
    data.frame(key = setdiff(a$key[a$SpD %in% spds], common), reason = rep("missing_target", length(setdiff(a$key[a$SpD %in% spds], common)))),
    data.frame(key = setdiff(b$key[b$SpD %in% spds], common), reason = rep("missing_source", length(setdiff(b$key[b$SpD %in% spds], common)))),
    data.frame(key = d$key[!valid], reason = rep("missing_covariate_or_too_few_spots", sum(!valid))))
  d <- d[valid, , drop = FALSE]
  rownames(d) <- d$key
  list(samples = d, exclusions = exclusions, design = make_design(d, cfg))
}
prepare_dge <- function(object, keys, design, restrict_genes = NULL) {
  counts <- object$counts[, match(keys, colnames(object$counts)), drop = FALSE]
  assert(identical(colnames(counts), keys), "Count alignment failed")
  if (!is.null(restrict_genes)) counts <- counts[rownames(counts) %in% restrict_genes, , drop = FALSE]
  dge <- edgeR::DGEList(counts)
  keep <- edgeR::filterByExpr(dge, design = design)
  assert(sum(keep) > 1L, "No usable gene universe")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  edgeR::calcNormFactors(dge, method = "TMM")
}
fit_voom <- function(dge, design, samples) {
  assert(identical(colnames(dge), samples$key) && identical(rownames(design), samples$key), "Model sample alignment failed")
  warns <- character()
  fit <- withCallingHandlers(edgeR::voomLmFit(dge, design = design, block = samples$donor,
    adaptive.span = TRUE, sample.weights = TRUE, keep.EList = TRUE),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  assert(!any(grepl("not estimable|setting to zero|Partial NA", warns, ignore.case = TRUE)), paste("Unstable fit:", paste(warns, collapse = "; ")))
  attr(fit, "screen_warnings") <- warns
  fit
}
cached_voom <- function(dge, design, samples, env) {
  key <- hash(list(env$signature, dge, design, samples$donor, "voom_adaptive_sampleweights"))
  path <- file.path(env$out, "cache", paste0(key, ".rds"))
  if (file.exists(path)) {
    entry <- readRDS(path)
    assert(identical(entry$key, key), "Corrupt model cache")
    fit <- entry$fit
  } else {
    fit <- fit_voom(dge, design, samples)
    save_atomic(list(key = key, fit = fit), path)
  }
  attr(fit, "cache_key") <- key
  fit
}
coef_table <- function(fit, term = "diagnosis") {
  contrast <- setNames(rep(0, ncol(fit$coefficients)), colnames(fit$coefficients))
  if (term == "diagnosis") {
    assert(all(c("DxSCZ", "DxNTC") %in% names(contrast)), "Diagnosis contrast unavailable")
    contrast[c("DxSCZ", "DxNTC")] <- c(1, -1)
  } else {
    assert(term %in% names(contrast), paste("Coefficient unavailable:", term))
    contrast[term] <- 1
  }
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts = matrix(contrast, ncol = 1)))
  tab <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none", confint = TRUE)
  assert(all(is.finite(tab$logFC)) && all(is.finite(tab$P.Value)), "Nonestimable gene coefficients")
  data.frame(gene_id = rownames(tab), beta = tab$logFC,
    se = as.numeric(fit$stdev.unscaled[, 1] * sqrt(fit$s2.post)),
    lower = tab$CI.L, upper = tab$CI.R, t = tab$t, p = tab$P.Value, q = tab$adj.P.Val)
}
prefix_coef <- function(tab, prefix) {
  names(tab)[-1] <- paste0(prefix, "_", names(tab)[-1]); tab
}
classify_pairs <- function(r, cfg) {
  r$historical_gate <- !is.na(r$historical_p) & r$historical_p < cfg$dx_p
  r$baseline_gate <- r$c_p < cfg$dx_p
  r$mediator_gate <- r$a_p < cfg$mediator_dx_p
  r$mediator_outcome_gate <- r$b_q < cfg$mediator_q
  r$significance_loss <- r$cprime_p >= cfg$dx_p
  gates <- c("historical_gate", "baseline_gate", if (!isFALSE(cfg$require_mediator_gate)) "mediator_gate",
             "mediator_outcome_gate", "significance_loss")
  r$screen_hit <- Reduce(`&`, r[gates])
  r$absolute_shrinkage <- abs(r$c_beta) - abs(r$cprime_beta)
  r$direction_reversal <- r$c_beta * r$cprime_beta < 0
  r$coefficient_shrinkage <- r$absolute_shrinkage > 0 & !r$direction_reversal
  r$relative_shrinkage <- ifelse(abs(r$c_beta) > 1e-8, 1 - r$cprime_beta / r$c_beta, NA_real_)
  r$attenuated_still_significant <- r$historical_gate & r$baseline_gate & r$coefficient_shrinkage & !r$significance_loss
  r$higher_priority <- r$screen_hit & r$coefficient_shrinkage & !is.na(r$b_q_global) & r$b_q_global < cfg$mediator_q
  r$failed_gates <- apply(r[, gates, drop = FALSE], 1, function(z) paste(gates[!z], collapse = ";"))
  r
}
assemble_results <- function(screen, historical, genes, source_fit, baseline, joint, cfg) {
  a <- coef_table(source_fit); a <- a[match(screen$mediator_id, a$gene_id), , drop = FALSE]
  assert(nrow(a) == 1 && !is.na(a$gene_id), "Mediator unavailable after expression filtering")
  r <- merge(prefix_coef(coef_table(baseline), "c"), prefix_coef(coef_table(joint), "cprime"), by = "gene_id", sort = FALSE)
  r <- merge(r, prefix_coef(coef_table(joint, "M"), "b"), by = "gene_id", sort = FALSE)
  r$gene_name <- genes$gene_name[match(r$gene_id, genes$gene_id)]
  for (field in names(a)[-1]) r[[paste0("a_", field)]] <- a[[field]]
  r$historical_p <- historical$p_value_scz[match(r$gene_id, historical$gene_id)]
  r$historical_q <- historical$fdr_scz[match(r$gene_id, historical$gene_id)]
  r$historical_beta <- historical$logFC_scz[match(r$gene_id, historical$gene_id)]
  for (field in c("screen_id", "source", "target", "mediator_id", "mediator_symbol")) r[[field]] <- screen[[field]]
  r$a_direction_agrees <- sign(r$a_beta) == screen$expected_direction
  r$same_gene_across_contexts <- r$gene_id == screen$mediator_id
  r$b_q_global <- NA_real_
  classify_pairs(r, cfg)
}
prepare_logcounts <- function(object, keys, restrict_genes = NULL) {
  E <- object$logcounts[, match(keys, colnames(object$logcounts)), drop = FALSE]
  assert(identical(colnames(E), keys), "Logcounts alignment failed")
  if (!is.null(restrict_genes)) E <- E[rownames(E) %in% restrict_genes, , drop = FALSE]
  E <- E[apply(E, 1, function(x) diff(range(x)) > 0), , drop = FALSE]
  assert(nrow(E) > 1L, "No usable gene universe")
  E
}
## manuscript engine: stored logcounts, donor correlation from the covariate design without SpD
## (spatialLIBD::registration_block_cor), one target correlation shared by nested models.
run_screen_logcounts <- function(screen, src, tgt, matched, cfg, restrict, extra) {
  d <- matched$samples
  des <- if (length(extra)) make_design(d, cfg, extra) else matched$design
  cd <- make_design(d, cfg, extra, include_spd = FALSE)
  sE <- prepare_logcounts(src, d$key, restrict$source)
  tE <- prepare_logcounts(tgt, d$key, restrict$target)
  assert(screen$mediator_id %in% rownames(sE), "Mediator unavailable in source logcounts")
  rho_s <- limma::duplicateCorrelation(sE, cd, block = d$donor)$consensus.correlation
  rho_t <- limma::duplicateCorrelation(tE, cd, block = d$donor)$consensus.correlation
  assert(is.finite(rho_s) && is.finite(rho_t), "Nonestimable donor correlation")
  sf <- limma::lmFit(sE, des, block = d$donor, correlation = rho_s)
  med <- as.numeric(sE[screen$mediator_id, ])
  assert(all(is.finite(med)) && sd(med) > 1e-8, "Missing or constant mediator expression")
  scaling <- c(mean = mean(med), sd = sd(med))
  d$M_logCPM <- med; d$M <- (med - scaling[["mean"]]) / scaling[["sd"]]
  jd <- make_design(d, cfg, c(extra, "M"))
  bf <- limma::lmFit(tE, des, block = d$donor, correlation = rho_t)
  jf <- limma::lmFit(tE, jd, block = d$donor, correlation = rho_t)
  bf$EList <- list(E = tE, weights = NULL)
  list(d = d, des = des, sdge = sE, tdge = tE, sf = sf, bf = bf, jf = jf, scaling = scaling,
       rho = c(source = rho_s, target = rho_t))
}
## fit an additional target design with the engine and variance settings of a primary artifact.
fit_target <- function(a, design, env) {
  if (is.null(a$rho)) return(cached_voom(a$tdge, design, a$samples, env))
  limma::lmFit(a$tdge, design, block = a$samples$donor, correlation = a$rho[["target"]])
}
run_screen <- function(screen, inputs, cfg, env, spds = cfg$spds, restrict = NULL, extra = character()) {
  src <- inputs$objects[[screen$source]]; tgt <- inputs$objects[[screen$target]]
  matched <- match_samples(src, tgt, cfg, spds)
  if (identical(cfg$engine, "limma_logcounts")) {
    z <- run_screen_logcounts(screen, src, tgt, matched, cfg, restrict, extra)
    r <- assemble_results(screen, inputs$historical[[screen$target]], tgt$genes, z$sf, z$bf, z$jf, cfg)
    r$n_samples <- nrow(z$d); r$n_donors <- length(unique(z$d$donor)); r$run_signature <- env$signature
    return(list(screen = screen, samples = z$d, exclusions = matched$exclusions, design = z$des,
                sdge = z$sdge, tdge = z$tdge, source_fit = z$sf, baseline = z$bf, joint = z$jf,
                scaling = z$scaling, rho = z$rho, results = r, signature = env$signature))
  }
  assert(identical(cfg$engine, "voom") && !length(extra), "Unsupported engine or extra covariates for voom run_screen")
  d <- matched$samples; des <- matched$design
  sdge <- prepare_dge(src, d$key, des, if (!is.null(restrict)) restrict$source else NULL)
  tdge <- prepare_dge(tgt, d$key, des, if (!is.null(restrict)) restrict$target else NULL)
  assert(screen$mediator_id %in% rownames(sdge), "Mediator failed matched expression filter")
  sf <- cached_voom(sdge, des, d, env)
  med <- as.numeric(sf$EList$E[screen$mediator_id, ])
  assert(all(is.finite(med)) && sd(med) > 1e-8, "Missing or constant mediator expression")
  scaling <- c(mean = mean(med), sd = sd(med))
  d$M_logCPM <- med; d$M <- (med - scaling[["mean"]]) / scaling[["sd"]]
  jd <- make_design(d, cfg, "M")
  bf <- cached_voom(tdge, des, d, env); jf <- cached_voom(tdge, jd, d, env)
  r <- assemble_results(screen, inputs$historical[[screen$target]], tgt$genes, sf, bf, jf, cfg)
  r$n_samples <- nrow(d); r$n_donors <- length(unique(d$donor)); r$run_signature <- env$signature
  list(screen = screen, samples = d, exclusions = matched$exclusions, design = des,
       sdge = sdge, tdge = tdge, source_fit = sf, baseline = bf, joint = jf,
       scaling = scaling, results = r, signature = env$signature)
}
fit_logcounts <- function(E, design, samples, correlation_design = design, weights = NULL) {
  rho <- limma::duplicateCorrelation(E, correlation_design, block = samples$donor, weights = weights)$consensus.correlation
  assert(is.finite(rho), "Nonestimable donor correlation")
  list(fit = limma::lmFit(E, design, block = samples$donor, correlation = rho, weights = weights), rho = rho)
}
fixed_comparison <- function(a, cfg) {
  E <- a$baseline$EList
  rho <- limma::duplicateCorrelation(E$E, a$design, block = a$samples$donor, weights = E$weights)$consensus.correlation
  assert(is.finite(rho), "Nonestimable fixed correlation")
  b <- limma::lmFit(E$E, a$design, block = a$samples$donor, correlation = rho, weights = E$weights)
  j <- limma::lmFit(E$E, make_design(a$samples, cfg, "M"), block = a$samples$donor, correlation = rho, weights = E$weights)
  list(baseline = b, joint = j, rho = rho)
}
