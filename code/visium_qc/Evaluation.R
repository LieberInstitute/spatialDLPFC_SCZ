# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(scater)
# library(pryr)                 # Check spe size
library(spatialLIBD)
library(tidyverse)


path_raw_spe <- here("processed-data/rds/spe",
                     "01_build_spe/", "spe_raw.rds")

raw_spe <- readRDS(
  path_raw_spe
)

spe <- raw_spe
# TODO: remove this
spe <- logNormCounts(spe)
spe <- logNormCounts(spe)

# TODO: how many proportion of the spots have missing




# TODO: Create pca
spe <- runPCA(spe) # TODO: edit this
plotPCA(spe, color_by = "sample_id")

spe <- runUMAP(spe) # TODO: edit this
plotUMAP(spe, color_by = "sample_id") +
  scale_alpha()
