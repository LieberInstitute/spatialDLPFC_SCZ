# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ComplexHeatmap)
  library(sessioninfo)
  library(viridis)
})



# Load Dx_DEG data ---
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
)

sig_gene <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    "Test_90DEGs.xlsx"
  ),
  col_names = FALSE
)[[1]]

sig_gene_df <- gene_df |>
  filter(gene %in% sig_gene)

ann_df <- sig_gene_df |>
  column_to_rownames(var = "gene") |>
  transmute(
    SCZ_reg = factor(
      logFC_scz > 0,
      levels = c(TRUE, FALSE),
      labels = c("Up", "Down")
    )
  )

# Load PB Data ----
spd_pb <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spatialDLPFC_BayesSpace_K09_pb_data.rds"
  )
)

overlap_genes <- intersect(sig_gene_df$ensembl, rownames(spd_pb))

overlap_genes |> length() # 86

missing_genes <- setdiff(sig_gene_df$ensembl, rownames(spd_pb))

gene_mat <- logcounts(spd_pb)[
  overlap_genes,
]


gene_mat_long <- gene_mat |>
  data.frame() |>
  rownames_to_column(var = "gene") |>
  pivot_longer(
    cols = starts_with("Br"),
    names_to = c("sample_id", "pos", "spd"),
    names_pattern = "^(Br\\d{4}_(mid|ant|post))_(.*)"
  )

miss_genes_mat <- matrix(
  NA,
  nrow = length(missing_genes), ncol = length(unique(gene_mat_long$spd))
)
rownames(miss_genes_mat) <- missing_genes

gene_mat_median <- gene_mat_long |>
  group_by(gene, spd) |>
  summarize(median = median(value)) |>
  ungroup() |>
  pivot_wider(
    id_cols = "gene",
    names_from = "spd",
    values_from = "median"
  ) |>
  column_to_rownames("gene") |>
  data.matrix() |>
  rbind(miss_genes_mat)

# Process PB Data

# Convert Ensembl to Gene name
stopifnot(nrow(gene_mat_median) == nrow(sig_gene_df))

rownames(gene_mat_median) <- sig_gene_df$gene[
  match(
    rownames(gene_mat_median),
    sig_gene_df$ensembl
  )
]

# Format SPD names
spd_name_df <- read_csv(
  here(
    "code/analysis/visium_spatial_clustering/",
    "bayesSpace_layer_annotations.csv"
  )
)


colnames(gene_mat_median) <- spd_name_df$layer_combo2[match(
  colnames(gene_mat_median),
  spd_name_df$cluster
)]





gene_names_hc_ordered <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)


pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_pec_spd_median_logCPM.pdf"
  ),
  height = 20
)
ComplexHeatmap::pheatmap(
  mat = gene_mat_median[
    gene_names_hc_ordered,
    order(colnames(gene_mat_median))
  ],
  name = "Scaled median logCPM",
  color = viridis(100, option = "magma"),
  scale = "row",
  column_title = "PEC SpD",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  row_split = ann_df[gene_names_hc_ordered, ],
  annotation_row = ann_df[gene_names_hc_ordered, , drop = FALSE],
  show_row_dend = FALSE,
  cellwidth = 10,
  cellheight = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_colors = list(SCZ_reg = c("Up" = "red", "Down" = "blue"))
)
dev.off()


# Session Info ----
session_info()
