# Job Script Related Notes ------------------------------------------------
# Running time: 1 hour
# Mem usage: asked for 40G, used 32GB in interactive session.

# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(sessioninfo)
})


# TODO: load in spe object
raw_spe <- readRDS(
  here::here("processed-data", "rds", "01_build_spe",
             # TODO: rename
             "test_raw_spe_w_spg_N63.rds")
)

# spe <- raw_spe[, raw_spe$in_tissue == 1 ]



# Subsetting Test Dataset -------------------------------------------------
test_spe <- raw_spe[, raw_spe$sample_id %in%
                      paste0("V13F27-294_",LETTERS[1:4],1)]

# Make sure all spots in test data have spg data.
stopifnot(
  ncol(test_spe) == test_spe |> colData() |> data.frame() |> 
  select(starts_with("spg_")) |> drop_na() |> 
  nrow()
  )

test_spe <- test_spe[, test_spe$in_tissue == TRUE]

saveRDS(test_spe,
        file = here("processed-data/rds/01_build_spe/test_spe_w_spg_n4.rds"))



# Session Info ----
session_info()
