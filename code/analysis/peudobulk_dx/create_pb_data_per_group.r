# Load library ----
library(SpatialExperiment)
library(tidyverse)
library(here)
library(sessioninfo)

# Load data ----
pb_sce <- read_rds(
  here(
    "processed-data", "rds", "layer_spd",
    # Focus on PRECAST 07 results
    # TODO: rename
    "test_spe_pseudo_PRECAST_07.rds"
  )
)


# Subset and create rds file for each dx group ----
pb_sce$dx |>
  unique() |>
  walk(
    ~ pb_sce[, pb_sce$dx == .x] |>
      saveRDS(
        file = here(
          "processed-data", "rds", "layer_spd",
          paste0("test_spe_pseudo_PRECAST_07_", .x, ".rds")
        )
      )
  )

# Session Info ----
sessioninfo::session_info()
