# Load packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(SingleR)
# library(tidyverse)


# Load spatialDLPFC data ---------------------------------------------------

spe <- fetch_data(type = "spatialDLPFC_Visium")

colData(spe) |> names()

set.seed(1)
scr_idx <- sample(unique(spe$sample_id), 15)
tar_idx <- setdiff(unique(spe$sample_id), scr_idx)

src_spe <- spe[,spe$sample_id %in% scr_idx]
tar_spe <- spe[,spe$sample_id %in% tar_idx]

# Spatial Registration ----------------------------------------------------

## Get enrichment statistics from source ---------------------------------
sce_modeling_results <- registration_wrapper(
  sce = src_spe,
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

