# Job Script Related Notes ------------------------------------------------
# Running time: 1 hour
# Mem usage: asked for 40G, used 32GB in interactive session.



# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(tidyverse)
library(sessioninfo)
library(spatialLIBD)


# TODO: load in spe object
raw_spe <- readRDS(
  here::here("processed-data", "rds", "01_build_spe",
             # TODO: rename
             "test_raw_spe_w_spg_N63.rds")
)

spe <- raw_spe[, raw_spe$in_tissue == 1 ]



# Visualize the distribution per sample -----------------------------------


c("spg_NBW", "spg_PBW", "spg_CNBW" ) |> 
  lapply(
    FUN = function(.x)
      vis_grid_gene(
      spe,
      geneid = .x,
      spatial = FALSE,
      point_size = 0.8,
      pdf_file = here(
        "plots/pnn",
        paste0("test_", .x, ".pdf")
      )
    )
  )




# Session Info ------------------------------------------------------------
session_info()