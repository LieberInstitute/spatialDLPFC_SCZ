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
## For a given program's top-2000 gene loadings, compute a loading-weighted
## average of log-normalized expression per spot:
##   score_j = sum_g(w_g * x_gj) / sum_g(w_g)
## restricted to genes matched (by gene symbol) on the Visium panel.
project_snap_score <- function(spe, loading_df, min_genes = 50) {
  program_name <- unique(loading_df$program)
  stopifnot(length(program_name) == 1)

  gene_match <- match(loading_df$gene, rowData(spe)$gene_name)
  keep <- !is.na(gene_match)

  cat(
    sprintf(
      "%s: matched %d / %d top-2000 genes on the Visium panel\n",
      program_name, sum(keep), nrow(loading_df)
    )
  )

  stopifnot(
    "Too few matched genes to compute a reliable projection" =
      sum(keep) >= min_genes
  )

  w <- loading_df$loading[keep]
  expr_mat <- logcounts(spe)[gene_match[keep], , drop = FALSE]

  as.numeric((w %*% expr_mat) / sum(w))
}

top2000 <- snap_loadings$top2000

snapa_score <- project_snap_score(
  spe,
  top2000 |> filter(program == "SNAP-a")
)

snapn_score <- project_snap_score(
  spe,
  top2000 |> filter(program == "SNAP-n")
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
  SNAPa_score = snapa_score,
  SNAPn_score = snapn_score
) |>
  mutate(
    # Normalized (z-scored) projection, computed across all spots/donors so
    # that SNAP-a and SNAP-n are on a comparable scale for visualization.
    SNAPa_zscore = as.numeric(scale(SNAPa_score)),
    SNAPn_zscore = as.numeric(scale(SNAPn_score))
  )

# error prevention
stopifnot(nrow(score_df) == ncol(spe))
stopifnot(!anyNA(score_df$SNAPa_zscore))
stopifnot(!anyNA(score_df$SNAPn_zscore))


# Save -----------------------------------------------------------------------
saveRDS(
  score_df,
  file.path(fld_out, "SNAP_score_df.rds")
)

print("Finished projecting SNAP-a / SNAP-n scores onto Visium spots")


# Session info ----
sessioninfo::session_info()
