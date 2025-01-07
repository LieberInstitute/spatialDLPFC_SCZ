# module load conda_R/4.3.x
# Note: to use the same version of PRECAST


# Load Packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(PRECAST)
  library(sessioninfo)
})

# Load PRECAST data ----
PRECASTObj <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_PRECASTObj_semi_inform_k_2-16.rds"
  )
)

PRECASTObj_model_select <- SelectModel(PRECASTObj)


# NOTE:
# Boyi deoesn't really like this step, in the sense that the batch correction step is slightly iffy for the following reasons
# 1. threatically speaking, the batch correction is done with house keeping genes. 
# How often the author would update the house keeping genes and where to acquire it.

# Recommentation
# It would probably better to export the PRECAST latent embedding, and run batch correction methods external to this step. 
seuInt <- IntegrateSpaData(PRECASTObj_model_select, species = "Human")


# TODO: run SelectModel to generate that sucker matrix for integration step

seuInt <- AddUMAP(seuInt,
  seed = 1
)
seuInt <- AddTSNE(seuInt, seed = 1)

# Calculate TSNE plot ----


#









## Example Integrated data ----



# seuInt <- readRDS(
#   here(
#     "processed-data/rds/spatial_cluster/PRECAST_prelim",
#     "test_seuIntObj_semi_inform_K8.rds"
#   )
# )



# seuInt
# Find dimensiona reduction names

dimPlot(seuInt,
  item = "cluster", reduction = "UMAP3",
  point_size = 0.0001, font_family = "serif"
)
dimPlot(seuInt,
  item = "batch", , reduction = "UMAP3",
  point_size = 0.0001, font_family = "serif"
)



dimPlot(seuInt,
  item = "cluster", reduction = "tSNE3",
  point_size = 0.0001, font_family = "serif"
)

dimPlot(seuInt,
  item = "batch", reduction = "tSNE3",
  point_size = 0.0001, font_family = "serif"
)


SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 1, combine = TRUE)


saveRDS(
  seuInt,
  here(
    "processed-data/rds/spatial_cluster/PRECAST",
    "test_seuInt_UMAP_tsne.rds"
  )
)


# Session Info ----
session_info()
