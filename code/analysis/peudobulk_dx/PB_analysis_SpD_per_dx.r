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

spd_rds <- list.files(
  here(
    "processed-data", "rds", "layer_spd"
  ),
  pattern = "(ntc|scz).rds$"
)


# Layer Enrichment analysis ----
for (.file in spd_rds) {
  # .file <- spd_rds[1]

  .spd <- str_remove(.file, "test_spe_pseudo_") |>
    str_remove("_(ntc|scz).rds")

  .dx_group <- str_remove(.file, "test_spe_pseudo_") |>
    str_remove(paste0(.spd, "_")) |>
    str_remove(".rds")

  sce_pseudo <- readRDS(
    here(
      "processed-data", "rds", "layer_spd",
      .file
    )
  )
  ## Calculate enrichment statistics ----
  layer_mod <-
    registration_model(
      sce_pseudo,
      covars = c("age", "sex", "slide_id"),
      var_registration = .spd
    )

  layer_block_cor <-
    registration_block_cor(
      sce_pseudo,
      registration_model = layer_mod,
      var_sample_id = "sample_id"
    )

  layer_res <- registration_stats_enrichment(
    sce_pseudo,
    block_cor = layer_block_cor,
    covars = c("age", "sex", "slide_id"),
    var_registration = .spd,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )

  # TODO: save analysis results
  saveRDS(
    layer_res,
    here(
      "processed-data", "rds", "layer_enrich_test",
      paste0("test_enrich_", .spd, "_", .dx_group, ".rds")
    )
  )
}

# Session Info ----
sessioninfo::session_info()
