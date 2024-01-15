# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(sessioninfo)
})


# Path --------------------------------------------------------------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)


# Load spe ----------------------------------------------------------------

spe <- readRDS(path_PRECAST_int_spe)


paste0("PRECAST_", seq.int(2,16)) |> 
walk(
  .f = function(.x){
    cat("Start for ", .x, "\n")
    vis_grid_clus(
    spe,
    clustervar = .x,
    point_size = 0.8,
    spatial = FALSE,
    alpha = 1,
    pdf_file = here("plots/spatial_cluster",
    paste0("test_semi_supervised_", .x, ".pdf"))
  )
    cat(.x, " finished\n")
  }
)




# SessionInfo -------------------------------------------------------------
sessioninfo::session_info()