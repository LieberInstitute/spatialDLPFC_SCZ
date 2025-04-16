# Objective:
# Check if the NAN value of the proportion and mean intensity. ISPG variable in the dataset to ensure it is correctly populated and analyze its distribution.

# Load lbrary ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
  library(here)
  library(sessioninfo)
})

# Load SPE object ----
spe <- readRDS(
  here::here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

col_df <- colData(spe) |> data.frame()


# Data Analysis ----
## DAPI channel ----
col_df$spg_PDAPI |> summary()

col_df |> filter(is.na(spg_IDAPI)) |>
pull(spg_PDAPI) |> summary()

col_df$spg_IDAPI |> summary()



## WFA channel ----
col_df$spg_PWFA |> summary()

# Session info ----
session_info()
