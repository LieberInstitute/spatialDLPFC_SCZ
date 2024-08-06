suppressPackageStartupMessages({
  library("DeconvoBuddies")
  library("SummarizedExperiment")
  library("SingleCellExperiment")
  library(readxl)
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load Pseudobulked data ----
spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

# Run Deconvo Buddies -----
## Deconvo buddies -----
tmp <- get_mean_ratio(
  spe_pb,
  cellType_col = "PRECAST_07",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)



# Find Marker genes
gene_df <- tmp |>
  group_by(cellType.target) |>
  slice_head(n = 5) |> ungroup()

gene_names <- gene_df$gene_ensembl

# Save genes
write_csv(
  tmp,
  here(
    "processed-data",
    "test_xenium_design_layer_genes_deconvo_buddies_all.csv"
  )
)


# Visualization -----
sub_samples <- read_xlsx(
  here::here("code/xenium_panel_design/Xenium_DonorList.xlsx"),
  col_names = FALSE
)[, 1:2] |> unlist()

sce <- spe_pb[gene_names, spe_pb$BrNumbr %in% sub_samples]


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
  annotation_row = gene_df |> column_to_rownames("gene_ensembl") |> select(cellType.target)
  # annotation_row = gene_df_raw |>
  #   filter(spatial_domain_resolution == "Sp09") |>
  #   group_by(test) |>
  #   arrange(fdr, .by_group = TRUE) |>
  #   slice_head(n = n_marker_gene) |>
  #   select(test, ensembl) |> rownames_to_column("ensembl")
)

# Session Info ----
session_info()
