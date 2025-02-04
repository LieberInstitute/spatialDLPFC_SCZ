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

# Save CSV file for Sang Ho ----
.spd <- "PRECAST_07"
layer_res <- readRDS(
  here(
    # TODO: need organize
    "processed-data", "rds", "layer_enrich_test",
    paste0("test_pairwise_", .spd, ".rds")
  )
)

# note: based on the distribution of p-value,
#  it feels this test is misspecificed.
hist(layer_res$`p_value_spd01-spd04`)

layer_res |>
  select(gene, ensembl, ends_with("spd01-spd04")) |>
  filter(`fdr_spd01-spd04` <= 0.05) |>
  arrange(`fdr_spd01-spd04`) |>
  tail()

## raw csv file ----
layer_res |> 
  write_csv(
    here(
      "code/analysis/04_SpD_marker_genes",
      "pairwise_DEG_enriched_all.csv"
    )
  )

## Spd01 genes, pos t_stat ----
layer_res |>
  select(gene, ensembl, ends_with("spd01-spd04")) |>
  filter(`fdr_spd01-spd04` <= 0.05, `t_stat_spd01-spd04` > 0) |>
  arrange(desc(`t_stat_spd01-spd04`)) |>
  write_csv(
    here(
      "code/analysis/04_SpD_marker_genes",
      "pairwise_DEG_enriched_in_spd01.csv"
    )
  )

layer_res |>
  select(gene, ensembl, ends_with("spd01-spd04")) |>
  filter(`fdr_spd01-spd04` <= 0.05, `t_stat_spd01-spd04` < 0) |>
  arrange(`t_stat_spd01-spd04`) |>
  write_csv(
    here(
      "code/analysis/04_SpD_marker_genes",
      "pairwise_DEG_enriched_in_spd04.csv"
    )
  )


# Session Info ----
session_info()
