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


# Load/Find() HVG ---------------------------------------------------------
# # Create a list of sample-specific SPEs
# spe_lst <- unique(spe$sample_id) |> 
#   set_names() |> 
#   map(~spe[,spe$sample_id == .x])
# 
# # Feed the list to multiBatchNorm
# tmp <- multiBatchNorm(spe_lst)
# 
# ## What is the diff between `multiBatchNorm` and lognorm
# # TODO
# 
# 
# # Create modelGeneVar for each item of the list
# 
# HVG_lst <- spe_lst |> map(.f = modelGeneVar)
# 
# # Use combineVar to combine these variables
# combined.dec <- combineVar(HVG_lst)


fld_ftr_slct <- here("processed-data/rds/feature_selection")
dec <- readRDS(here(fld_ftr_slct, "HVG_imp2.rds"))

chosen.hvgs <- getTopHVGs(dec, n=5000)

# Validation of batch effect
combined <- correctExperiments(spe_lst,
                               PARAM=NoCorrectParam()) # Note not corrected 


# set.seed(100)
# combined <- runPCA(combined, subset_row=chosen.hvgs)
# combined <- runTSNE(combined, dimred="PCA")
# before_plot <- plotTSNE(combined, colour_by="batch")
# ggplot2::ggsave(here::here("plots/test_tsne_before_MNN_correct.pdf"), before_plot)



# Fast MNN ----------------------------------------------------------------
set.seed(101)
# NOTE: is there a diff if I use spe, instead of combine from `correctExperiments`
f.out <- fastMNN(combined, batch=combined$batch, subset.row=chosen.hvgs)

fld_batch_correct <- here("processed-data", "rds", "spe", "batch_corrected")
saveRDS(object = f.out, here(fld_batch_correct, "test_spe_MNN.rds") )


# set.seed(103)
# f.out <- runTSNE(f.out, dimred="corrected")
# after_plot <- plotTSNE(f.out, colour_by="batch")
# ggplot2::ggsave(here::here("plots/test_tsne_after_MNN_correct.pdf"), after_plot)

