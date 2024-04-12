# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(tidyverse)
  library(SpotSweeper)
  library(sessioninfo)
  library(here)
})

stopifnot(
  "To replicate the result, please use spotsweeper 0.99.2" =
    package.version("SpotSweeper") >= "0.99.2"
)

# Load Data ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

## Remove out tissue spots
ret_spe <- spe[, spe$in_tissue == TRUE]

# Spatially-aware QC (SpotSweeper) ----
# total umi
ret_spe <- localOutliers(
  ret_spe,
  metric = "sum_umi",
  direction = "lower",
  log = TRUE
)

# unique genes
ret_spe <- localOutliers(ret_spe,
  metric = "sum_gene",
  direction = "lower",
  log = TRUE
)

# mitochondrial percent
ret_spe <- localOutliers(
  ret_spe,
  metric = "expr_chrM_ratio",
  direction = "higher",
  log = FALSE
)

# Summarize metrics
ret_spe$local_outliers <- as.logical(ret_spe$sum_umi_outliers) |
  as.logical(ret_spe$sum_gene_outliers) |
  as.logical(ret_spe$expr_chrM_ratio_outliers)


# Create & save outlier df ----
outlier_df <- data.frame(
  key = spe$key,
  local_outliers = (spe$key %in% ret_spe$key[ret_spe$local_outliers == TRUE])
)

saveRDS(
  outlier_df,
  file = here::here(
    "processed-data/02_visium_qc",
     "local_outlier_df.rds"
  )
)

# Session info ----
session_info()
