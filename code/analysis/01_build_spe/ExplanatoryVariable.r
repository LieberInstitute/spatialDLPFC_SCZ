# SLurm info
# Time: 4hr
# Space: 65 GB

# Load Library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scater)
  library(sessioninfo)
})

# Load Data ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63.rds"
  )
)

# Pre-QC filtering ----
## Remove out tissue spots ----
spe <- spe[, spe$in_tissue == TRUE]
ncol(spe)
## Remove sum_umi=0 spots ----
# Remove these spots to do lognormalization
spe <- spe[, spe$sum_umi != 0]
ncol(spe)


# Preprocessing ----
## Fetch demo info
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$age <- metadata(spe)$dx_df$age[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

## Log transformation -----
# Create logcounts
spe <- logNormCounts(spe)

# Gene Variance Explained ----
## Note: not run well

var_mat <- getVarianceExplained(
  spe,
  variables = c("dx", "sex", "age", "sample_id")
)

saveRDS(var_mat,
        here("varExplained.rds"))

plotExplanatoryVariables(var_mat)

# Session ----
sessioninfo::session_info()