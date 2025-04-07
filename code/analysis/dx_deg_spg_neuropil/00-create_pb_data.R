# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  # library(limma)
  library(sessioninfo)
  library(here)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(
  here::here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

colData(spe) |> names()
