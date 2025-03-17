# Load Libray --------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(sessioninfo)
  library(tidyverse)
  library(ggrepel)
})

# Load data -----
# TODO: find the pseudobulk data

# Pseudo-bulk DE ----
spe_bulk <-
  registration_pseudobulk(
    spe,
    var_registration = "dx",
    var_sample_id = "sample_id",
    min_ncells = 10
  )

covars <- c("age", "sex", "lot_num")

set.seed(20240411)
spe_bulk <- runPCA(sce_pseudo)

tmp <- registration_stats_enrichment(
  spe_bulk,
  block_cor = NaN,
  covars = c("age", "sex"),
  # var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

tmp |>
  filter(fdr_scz <= 0.05) |>
  nrow()

tmp <- registration_stats_enrichment(
  spe_bulk,
  block_cor = NaN,
  covars = c("age", "sex", "lot_num"),
  # var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)



plotPCA(
  spe_bulk,
  colour_by = "spd",
  ncomponents = 6,
  point_size = 0.5,
  label_format = c("%s %02i", " (%i%%)"),
)





registration_mod <-
  registration_model(spe_bulk,
    covars = covars
  )


block_cor <-
  registration_block_cor(
    spe_bulk,
    registration_model = registration_mod
  )

results_enrichment <-
  registration_stats_enrichment(
    spe_bulk,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )

saveRDS(
  results_enrichment,
  file = here("processed-data/rds/pseudo_bulk/sample_bulk_DE.rds")
)


# Session Info ----
sessioninfo::session_info()
