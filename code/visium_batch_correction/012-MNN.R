library(batchelor)
library(scuttle)
library(here)
library(SpatialExperiment)
library(tidyverse)
library(scran)
library(scater)


# Read in spe
spe <- readRDS(here::here("processed-data", "rds", "spe",
                          "spe_after_spot_qc.rds"))


# Create a list of sample-specific SPEs
spe_lst <- unique(spe$sample_id) |> 
  set_names() |> 
  map(~spe[,spe$sample_id == .x])

# Feed the list to multiBatchNorm
tmp <- multiBatchNorm(spe_lst)

## What is the diff between `multiBatchNorm` and lognorm
# TODO


# Create modelGeneVar for each item of the list

HVG_lst <- spe_lst |> map(.f = modelGeneVar)

# Use combineVar to combine these variables
combined.dec <- combineVar(HVG_lst)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)

# Validation of batch effect
combined <- correctExperiments(spe_lst,
                               PARAM=NoCorrectParam()) # Note not corrected 


set.seed(100)
combined <- runPCA(combined, subset_row=chosen.hvgs)
combined <- runTSNE(combined, dimred="PCA")
plotTSNE(combined, colour_by="batch")


# Fast MNN
set.seed(101)
f.out <- fastMNN(combined, batch=combined$batch, subset.row=chosen.hvgs)


set.seed(103)
f.out <- runTSNE(f.out, dimred="corrected")
plotTSNE(f.out, colour_by="batch")