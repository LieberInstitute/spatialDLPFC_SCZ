# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(sessioninfo)
  library(tidyverse)
})



# Mem requested
# TODO: adjust
# First try 20G


# Config -----------------------------------------------------------
## Path Config -------------------------------------------------------------
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



# Pseudo-bulk DE ----------------------------------------------------------
## Add dx informaiton to the dataset
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]


## Subset Layers -----------------------------------------------------------

spe$spd <- spe$PRECAST_5

spe$spd <- spe$spd |> 
  factor(labels = paste0("SpD_", unique(spe$PRECAST_5)))

spd_list <- spe$spd |>
  levels() |> 
  map(.f = function(spe, .spd){
    spe[,spe$spd == .spd]
  },
  spe = spe)

## Run DE ---------------

spe_sub <- spd_list[[1]]
sce_pseudo <-
  registration_pseudobulk(
    spe_sub,
    var_registration = "dx",
    var_sample_id = "sample_id",
    min_ncells = 10
  )

# WHAT TODO with this?



# Session info ------------------------------------------------------------
sessioninfo::session_info()

