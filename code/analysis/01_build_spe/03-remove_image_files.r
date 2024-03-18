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

imgData(spe) <- NULL

saveRDS(
  spe,
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

# Report Session Info ----
sessioninfo::session_info()