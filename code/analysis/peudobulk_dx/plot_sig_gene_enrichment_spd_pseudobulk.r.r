# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SpatialExperiment)
  library(here)
  library(sessioninfo)
})

# Load Data ----

spe_pb <- readRDS(
  here(
    "processed-data/rds/layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

