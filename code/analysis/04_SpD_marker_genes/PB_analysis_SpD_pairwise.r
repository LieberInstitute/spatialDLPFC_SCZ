# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

# Load data -----
## SpD07 ----
.spd <- "PRECAST_07"

sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)


# Pariwise test ----
layer_mod <-
  registration_model(
    sce_pseudo,
    covars = c(
      "age", "sex" # , "slide_id"),
    ),
    var_registration = .spd
  )

layer_block_cor <-
  registration_block_cor(
    sce_pseudo,
    registration_model = layer_mod,
    var_sample_id = "sample_id"
  )

layer_res <- registration_stats_pairwise(
  sce_pseudo,
  block_cor = layer_block_cor,
  registration_model = layer_mod,
  # Comment out first
  # NOTE: not sure if this is correct to do statisitcally
  # covars = c("age", "sex", "slide_id"),
  var_registration = .spd,
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

layer_res |>
  arrange(`fdr_spd01-spd04`) |>
  select(`fdr_spd01-spd04`, `logFC_spd01-spd04`, gene) |>
  head(n = 20)

# Save pairwise test result ----
saveRDS(
  layer_res,
  here(
    # TODO: need organize
    "processed-data", "rds", "layer_enrich_test",
    paste0("test_pairwise_", .spd, ".rds")
  )
)


# Session Info ----
session_info()
