# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  # library(ggrepel)
  library(sessioninfo)
})

# Load donor-spd PB data ----
sce_pseudo <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_PRECAST07_donor_spd.rds"
  )
)

# Dx DE analysis ----
dx_mod <- registration_model(
  sce_pseudo,
  covars = c("fnl_spd", "age", "sex", "rin", "lot_num"),
  var_registration = "dx"
)

dx_block_cor <- registration_block_cor(
  sce_pseudo,
  registration_model = dx_mod,
  var_sample_id = "sample_id"
)

dx_res <- registration_stats_enrichment(
  sce_pseudo,
  block_cor = dx_block_cor,
  covars = c("fnl_spd", "age", "sex", "rin", "lot_num"),
  var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

# Save Outcome ----
## Save as RDS ----

## Save as CSV ----
dx_res |>
  arrange(fdr_scz) |>
  write.csv(
    here(
      "processed-data/rds/10_dx_deg_adjust_spd",
      "dx-deg_rin_adj_PRECAST07.csv"
    )
  )

# Session Info ----
sessioninfo::session_info()
