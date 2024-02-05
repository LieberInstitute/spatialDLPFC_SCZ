# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  # library(spatialLIBD)
  library(sessioninfo)
  library(tidyverse)
  # library(ggrepel)
})



# Config -----------------------------------------------------------
## Path Config -------------------------------------------------------------
# TODO: decide if this needs to be revised 
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

# Load data ---------------------------------------------------------------
spe <- readRDS(
  path_PRECAST_int_spe
)



# NOTE: For testing only
# spe_backup <- spe
# subset_id <- metadata(spe)$dx_df |> group_by(dx) |>
#   slice_head(n=2) |> ungroup() |>
#   pull(sample_id)
# 
# spe <- spe_backup[ ,spe_backup$sample_id %in% subset_id]

## Add dx informaiton to the dataset
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]

spe$age <- metadata(spe)$dx_df$age[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]

spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]



# Subsetting one dx group -------------------------------------------------
spe <- spe[, spe$dx == unique(spe$dx)[1]]

cat("Total number of spots is ", ncol(spe), "for ", unique(spe$dx)[1], ".\n")
# NMF ---------------------------------------------------------------------
library(RcppML)
library(Matrix)
# library(CoGAPS)
# library(projectR)
# library(scuttle)
# library(spatialLIBD)
# library(escheR)
# library(here)
# library(tidyverse)
# library(ggforce)

# Code adapted from Cindy's code
# https://github.com/cindyfang70/xenium-sandbox/blob/main/code/NMF/01_exploratory_visium_nmf.R


## Config ------------------------------------------------------------------
k <- 2

spe <- logNormCounts(spe)
A <- logcounts(spe) # using logcounts because there are multiple datasets

model <- RcppML::nmf(A, k = k, seed=1237)


## Save results
saveRDS(
  model, 
  file = here(
    "processed-data/rds/NMF",
    paste0("test_per_dx_", unique(spe$dx)[1], ".rds")
  )
)



# Session Info ------------------------------------------------------------

sessioninfo::session_info()
