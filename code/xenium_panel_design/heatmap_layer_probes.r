# Load packages ----
# suppressPackageStartupMessages({
#   library(here)
#   library(tidyverse)
#   library(SingleCellExperiment)
#   library(readxl)
#   library(ComplexHeatmap)
#   library(viridis)
#   library(sessioninfo)
# })


# n_marker_gene <- 5

# gene_df_raw <- read.csv(
#   here(
#     # TODO: change path
#     "code/analysis/visium_spatial_clustering",
#     "TableS8_sig_genes_FDR5perc_enrichment.csv"
#   )
# )

# gene_df <- gene_df_raw |>
#   filter(spatial_domain_resolution == "Sp09") |>
#   group_by(test) |>
#   arrange(fdr, .by_group = TRUE) |>
#   slice_head(n = n_marker_gene) |>
#   ungroup() |>
#   transmute(ensembl, spd = test)

# gene_df <- gene_df[duplicated(gene_df$ensembl) == FALSE, ]

# gene_names <- gene_df$ensembl


# # gene_df <- gene_df_raw |>
# #   filter(spatial_domain_resolution == "Sp09") |>
# #   group_by(test) |>
# #   arrange(fdr, .by_group = TRUE) |>
# #   slice_head(n = n_marker_gene) |>
# #   filter(ensembl %in% gene_names) |>
# #   transmute(ensembl, spd = spatial_domain_resolution)


# # Load PB Data ----
# ## All data ---
# spe_pb <- readRDS(
#   here(
#     "processed-data", "rds", "layer_spd",
#     "test_spe_pseudo_PRECAST_07.rds"
#   )
# )

# ## Subset to experiment data ----
# sub_samples <- read_xlsx(
#   here::here("code/xenium_panel_design/Xenium_DonorList.xlsx"),
#   col_names = FALSE
# )[, 1:2] |> unlist()





create_heatmap <- function(spe_pb, gene_names) {
  sce <- spe_pb[gene_names, ]
  col_df <- colData(sce) |>
    data.frame() |>
    select(PRECAST_07, dx, sample_id) |>
    mutate(
      PRECAST_07 = factor(PRECAST_07,
        levels = c(
          "spd07", "spd06", "spd02",
          "spd05", "spd03", "spd01", "spd04"
        )
      )
    ) |>
    arrange(
      PRECAST_07, dx, sample_id
    )

  # col_df

  gene_mat <- logcounts(sce)[
    ,
    rownames(col_df)
  ]


  ComplexHeatmap::pheatmap(
    mat = gene_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "row",
    annotation_col = col_df |> select(PRECAST_07),
    annotation_row = gene_df |> column_to_rownames("ensembl") |> select(spd)
    # annotation_row = gene_df_raw |>
    #   filter(spatial_domain_resolution == "Sp09") |>
    #   group_by(test) |>
    #   arrange(fdr, .by_group = TRUE) |>
    #   slice_head(n = n_marker_gene) |>
    #   select(test, ensembl) |> rownames_to_column("ensembl")
  )
}


# Session Info ----
# session_info()
