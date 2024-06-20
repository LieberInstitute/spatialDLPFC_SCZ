# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(sessioninfo)
})

# Load data ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "test_raw_spe_UMAP_N63.rds"
  )
)

## Fetch demo info
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$brnum <- metadata(spe)$dx_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sample_id <- paste0(
  spe$brnum, "_", spe$dx
)


# Session Info ---
session_info()
