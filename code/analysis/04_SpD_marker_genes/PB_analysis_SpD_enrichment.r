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
  pattern = ".rds"
)


# Layer Enrichment analysis ----
for (.file in spd_rds) {
  # .file <- spd_rds[1]

  .spd <- str_remove(.file, "test_spe_pseudo_") |>
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
      # TODO: need organize
      "processed-data", "rds", "layer_enrich_test",
      paste0("test_enrich_", .spd, ".rds")
    )
  )

  layer_res <- readRDS(here(
      # TODO: need organize
      "processed-data", "rds", "layer_enrich_test",
      paste0("test_enrich_PRECAST_07.rds")
    ))
  ## Format enrichment test res ----
  t_stats <- layer_res[, grep("^t_stat_", colnames(layer_res))]
  colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))

  ## Registration to Manual Annotationn ----
  manual_modeling_results <- fetch_data(type = "modeling_results")

  ### top 100 genes only ---
  manual_cor_g100 <- layer_stat_cor(
    t_stats,
    manual_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
  )
  pdf(
    here(
      "plots", "cluster_anno",
      paste0("test_", .spd, "_manual_g100.pdf")
    )
  )
  layer_stat_cor_plot(
    manual_cor_g100,
    max = max(manual_cor_g100)
  )
  dev.off()

  ### all genes ---
  manual_cor_all <- layer_stat_cor(
    t_stats,
    manual_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = NULL
  )
  pdf(
    here(
      "plots", "cluster_anno",
      paste0("test_", .spd, "_manual_all.pdf")
    )
  )
  layer_stat_cor_plot(
    manual_cor_all,
    max = max(manual_cor_all)
  )
  dev.off()

  ## Registration to spatialDLPFC ----
  spatialDLPFC_modeling_results <- fetch_data(
    type = "spatialDLPFC_Visium_modeling_results"
  )

  ### top 100 genes only ---
  spatialDLPFC_cor_g100 <- layer_stat_cor(
    t_stats,
    spatialDLPFC_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = 100
  )
  pdf(
    here(
      "plots", "cluster_anno",
      paste0("test_", .spd, "_spatialDLPFC_g100.pdf")
    )
  )
  layer_stat_cor_plot(
    spatialDLPFC_cor_g100,
    max = max(spatialDLPFC_cor_g100)
  )
  dev.off()

  ### all genes ---
  spatialDLPFC_cor_all <- layer_stat_cor(
    t_stats,
    spatialDLPFC_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = NULL
  )
  pdf(
    here(
      "plots", "cluster_anno",
      paste0("test_", .spd, "_spatialDLPFC_all.pdf")
    )
  )
  layer_stat_cor_plot(
    spatialDLPFC_cor_all,
    max = max(spatialDLPFC_cor_all)
  )
  dev.off()
}

# Session Info ----
sessioninfo::session_info()
