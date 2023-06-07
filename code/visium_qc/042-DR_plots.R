# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(scater)
# library(pryr)                 # Check spe size
# library(spatialLIBD)
library(tidyverse)


path_spe_after_spot_qc <- here::here("processed-data", "rds", "spe", "spe_after_spot_qc.rds")


spe <- readRDS(
  path_spe_after_spot_qc
)

# Create logcounts
spe <- logNormCounts(spe)

spe$slide_num <- str_remove(spe$sample_id, "_[A-D][1]")

spe$dx <- metadata(spe)$dx[match(spe$sample_id, metadata(spe)$sample_name)]


# TODO: Create pca
spe <- runPCA(spe) # TODO: edit this
plotPCA(spe, color_by = "sample_id")

plotPCA(spe, color_by = "slide_num")


set.seed(1)
spe <- runUMAP(spe) # TODO: edit this




plotUMAP(spe, color_by = "sample_id") 
plotUMAP(spe[, spe$sample_id="V12F14-057_A1"], color_by = "sample_id")
plotUMAP(spe, color_by = "slide_num")
plotUMAP(spe, color_by = "dx")
