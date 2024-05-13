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



## Subset NTC sample only -----------------------------------------------------------

sub_spe <- spe[, spe$dx == "ntc"]

# TODO: examine every PRECAST K
sub_spe$spd <- sub_spe$PRECAST_5

sub_spe$spd <- sub_spe$spd |> 
  factor(labels = paste0("SpD_", unique(spe$PRECAST_5)))

# spd_list <- spe$spd |>
#   levels() |> 
#   map(.f = function(spe, .spd){
#     spe[,spe$spd == .spd]
#   },
#   spe = spe)

## Create Pseudobulk data --------------
# # TODO: make the following code to iterate over the number of spds
# spe_sub <- spd_list[[2]]

# THis is taking the sum. The logcounts is CPM scale.
sce_pseudo <-
  registration_pseudobulk(
    sub_spe,
    var_registration = "spd",
    var_sample_id = "sample_id",
    min_ncells = 10
  )

covars <- c("age", "sex")
registration_mod <-
  registration_model(sce_pseudo,
                     covars = covars)

block_cor <-
  registration_block_cor(
    sce_pseudo,
    registration_model = registration_mod
  )

results_enrichment <-
  registration_stats_enrichment(
    sce_pseudo,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )

suffix <- "all"
results_anova <-
  registration_stats_anova(
    sce_pseudo,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    suffix = suffix
  )

results_enrichment_noadj <-
  registration_stats_enrichment(
    sce_pseudo,
    block_cor = block_cor,
    covars = NULL,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )

# How is the data summarized?


# https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/07_layer_differential_expression/03_model_BayesSpace.R

# WHAT TODO with this?



# Session info ------------------------------------------------------------
sessioninfo::session_info()

