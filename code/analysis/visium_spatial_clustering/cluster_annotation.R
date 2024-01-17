# Load packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
# library(tidyverse)


# Path --------------------------------------------------------------------
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

col_df <- colData(spe) |> data.frame()


# -------------------------------------------------------------------------
vars <- grep("^PRECAST", names(col_df), value = TRUE)

# TODO: write a loop
.var <- vars[[4]]

spe$spd <- paste0("SpD_", spe[[.var]])

# TODO: consider to ajdust for age, sex.
PNN_modeling_results <- registration_wrapper(
  sce = spe,
  var_registration = "spd",
  var_sample_id = "sample_id",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

str(PNN_modeling_results)

# Extract Enrichment t-statistics ------
PNN_t_stats <- PNN_modeling_results$enrichment[, grep("^t_stat", colnames(PNN_modeling_results$enrichment))]
colnames(PNN_t_stats) <- gsub("^t_stat_", "", colnames(PNN_t_stats))


# Down load spatialDLPFC modeling result ----------------------------------

DLPFC_modeling_results <- fetch_data("spatialDLPFC_Visium_modeling_results")

str(DLPFC_modeling_results)

DLPFC_t_stats <- DLPFC_modeling_results$enrichment[, grep("^t_stat", colnames(DLPFC_modeling_results$enrichment))]
colnames(DLPFC_t_stats) <- gsub("^t_stat_", "", colnames(DLPFC_t_stats))

# cor_layer <- layer_stat_cor(
cor_layer <- layer_stat_cor(
  stats = PNN_t_stats,
  modeling_results = DLPFC_modeling_results,
  model_type = "enrichment",
  top_n = 100
)

layer_stat_cor_plot(cor_layer, max = max(cor_layer))

