# Load packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(SingleR)
# library(tidyverse)


# Load spatialDLPFC data ---------------------------------------------------

spe <- fetch_data(type = "spatialDLPFC_Visium")

colData(spe) |> names()

# Error
# spe$BayesSpace_harmony_09 <- factor(spe$BayesSpace_harmony_09)
spe$BayesSpace_harmony_09 <- factor(
  spe$BayesSpace_harmony_09,
  labels = paste0("k09d", levels(spe$BayesSpace_harmony_09))
)
levels(spe$BayesSpace_harmony_09 )

set.seed(1)
scr_idx <- sample(unique(spe$sample_id), 5)
tar_idx <- setdiff(unique(spe$sample_id), scr_idx)

src_spe <- spe[,spe$sample_id %in% scr_idx]
tar_spe <- spe[,spe$sample_id %in% tar_idx]

# Spatial Registration ----------------------------------------------------

## Get enrichment statistics from source ---------------------------------
src_modeling_results <- registration_wrapper(
  sce = src_spe,
  # sce = src_spe,
  var_registration = "BayesSpace_harmony_09",
  var_sample_id = "sample_id",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

tar_modeling_results <- registration_wrapper(
  sce = tar_spe,
  # sce = src_spe,
  var_registration = "BayesSpace_harmony_09",
  var_sample_id = "sample_id",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)


# Error in combn(colnames(registration_model)[regis_cols], 2) : n < m

## Extract enrichment t-statistics  ---------------------------------


# Single R ----------------------------------------------------------------
predictions <- SingleR(test=tar_spe,
                       ref=src_spe, labels=src_spe$BayesSpace_harmony_09)

