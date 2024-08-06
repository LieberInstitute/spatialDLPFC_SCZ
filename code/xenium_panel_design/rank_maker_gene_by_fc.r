suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SpatialExperiment)
  library(sessioninfo)
  library(readxl)
})


n_marker_gene <- 5

gene_df_raw <- read.csv(
  here(
    # TODO: change path
    "code/analysis/visium_spatial_clustering",
    "TableS8_sig_genes_FDR5perc_enrichment.csv"
  )
)

unique()


setdiff(
  gene_df_raw |>
    filter(spatial_domain_resolution == "Sp09") |>
    group_by(test) |>
    arrange(fdr, .by_group = TRUE) |>
    slice_head(n = n_marker_gene) |>
    pull(gene),
  gene_df_raw |>
    filter(spatial_domain_resolution == "Sp09") |>
    group_by(test) |>
    arrange(desc(stat), .by_group = TRUE) |>
    slice_head(n = n_marker_gene) |>
    pull(gene)
)


setdiff(
  gene_df_raw |>
    filter(spatial_domain_resolution == "Sp09") |>
    group_by(test) |>
    arrange(desc(stat), .by_group = TRUE) |>
    slice_head(n = n_marker_gene) |>
    pull(gene),
  gene_df_raw |>
    filter(spatial_domain_resolution == "Sp09") |>
    group_by(test) |>
    arrange(fdr, .by_group = TRUE) |>
    slice_head(n = n_marker_gene) |>
    pull(gene)
)
