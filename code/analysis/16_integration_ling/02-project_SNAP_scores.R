# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(scater)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
fld_out <- here("processed-data/rds/16_integration_ling")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)


# Load data ----
## Finalized spe (63 donors, kept/QC-passed spots only) ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

## Ling et al. SNAP gene loadings, prepared in 01-prepare_SNAP_loadings.R ----
snap_loadings <- readRDS(
  file.path(fld_out, "SNAP_loadings.rds")
)


# Project SNAP loadings onto spots -----------------------------------------
## For a given program's gene loadings, compute a per-spot projection from
## log-normalized expression, restricted to genes matched (by gene symbol) on
## the Visium panel:
##   - top-2000 gene sets have all-nonnegative loadings (Ling et al.'s ranked-
##     by-loading convention), so we compute a loading-weighted *average*:
##       score_j = sum_g(w_g * x_gj) / sum_g(w_g)
##     Dividing by the total loading prevents scores from differing only
##     because one gene set has larger aggregate loading.
##   - full gene sets include negative loadings (genes pushed down by the
##     program, not just up), so a weighted-average normalization by sum(w)
##     is not meaningful (the denominator can be small, zero, or sign-
##     flipping). For these, we instead compute the raw loading-weighted sum:
##       score_j = sum_g(w_g * x_gj)
project_snap_score <- function(spe, loading_df, gene_set, min_genes = 50) {
  program_name <- unique(loading_df$program)
  stopifnot(length(program_name) == 1)

  # Safeguard against duplicate symbols in future loading-file releases.
  loading_df <- loading_df |>
    group_by(gene) |>
    summarize(loading = sum(loading), .groups = "drop")

  gene_match <- match(loading_df$gene, rowData(spe)$gene_name)
  keep <- !is.na(gene_match)

  cat(
    sprintf(
      "%s (%s): matched %d / %d genes on the Visium panel\n",
      program_name, gene_set, sum(keep), nrow(loading_df)
    )
  )

  stopifnot(
    "Too few matched genes to compute a reliable projection" =
      sum(keep) >= min_genes
  )

  w <- loading_df$loading[keep]
  expr_mat <- logcounts(spe)[gene_match[keep], , drop = FALSE]

  # Only normalize into a weighted average when all matched loadings are
  # nonnegative (as in the top-2000 gene sets); otherwise return the raw
  # loading-weighted sum, since normalizing by sum(w) is not meaningful when
  # loadings can be negative.
  if (all(w >= 0)) {
    stopifnot(
      "Matched gene loadings must have a nonzero sum" = sum(w) != 0
    )
    as.numeric((w %*% expr_mat) / sum(w))
  } else {
    as.numeric(w %*% expr_mat)
  }
}

top2000 <- snap_loadings$top2000
all_genes <- snap_loadings$full

snapa_score <- project_snap_score(
  spe,
  top2000 |> filter(program == "SNAP-a"),
  gene_set = "top 2000"
)

snapn_score <- project_snap_score(
  spe,
  top2000 |> filter(program == "SNAP-n"),
  gene_set = "top 2000"
)

# NOTE: the full loading lists contain negative loadings, so
# `project_snap_score()` returns a raw loading-weighted sum here (not
# normalized by sum(w), unlike the top-2000 scores above).
snapa_all_genes_score <- project_snap_score(
  spe,
  all_genes |> filter(program == "SNAP-a"),
  gene_set = "all genes"
)

snapn_all_genes_score <- project_snap_score(
  spe,
  all_genes |> filter(program == "SNAP-n"),
  gene_set = "all genes"
)


# Build slim, key-indexed score data frame ----------------------------------
# NOTE: kept slim (not re-saving the whole `spe`) so it can be `left_join`-ed
# back into `colData(spe)` in downstream scripts, following the convention
# used for PRECAST labels in `code/analysis/01_build_spe/20_create_final_spe.R`
score_df <- tibble(
  key = spe$key,
  sample_id = spe$sample_id,
  brnum = spe$brnum,
  dx = spe$dx,
  spd_label = spe$spd_label,
  SNAPa_top_2000_score = snapa_score,
  SNAPn_top_2000_score = snapn_score,
  SNAPa_all_genes_score = snapa_all_genes_score,
  SNAPn_all_genes_score = snapn_all_genes_score
)

# error prevention
stopifnot(nrow(score_df) == ncol(spe))
stopifnot(!anyNA(score_df$SNAPa_top_2000_score))
stopifnot(!anyNA(score_df$SNAPn_top_2000_score))
stopifnot(!anyNA(score_df$SNAPa_all_genes_score))
stopifnot(!anyNA(score_df$SNAPn_all_genes_score))


# Save -----------------------------------------------------------------------
saveRDS(
  score_df,
  file.path(fld_out, "SNAP_score_df.rds")
)

print("Finished projecting top-2000 and all-gene SNAP-a / SNAP-n scores onto Visium spots")


# Session info ----
sessioninfo::session_info()
