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
snRNA_pb <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spatialDLPFC_sn_hc_pb_data.rds"
  )
)

gene_mat <- logcounts(snRNA_pb)[sig_gene_df$ensembl, ]


gene_mat_long <- gene_mat |>
  data.frame() |>
  rownames_to_column(var = "gene") |>
  pivot_longer(
    cols = starts_with("Br"),
    names_to = c("sample_id", "pos", "cell_type"),
    names_pattern = "^(Br\\d{4}_(mid|ant|post))_(.*)"
  )

gene_mat_median <- gene_mat_long |>
  group_by(gene, cell_type) |>
  summarize(median = median(value)) |>
  ungroup() |>
  pivot_wider(
    id_cols = "gene",
    names_from = "cell_type",
    values_from = "median"
  ) |>
  column_to_rownames("gene") |>
  data.matrix()

# Process PB Data

## Function ---



# Need to accumulate genes




# Convert Ensembl to Gene name
rownames(gene_mat_median) <- rowData(snRNA_pb)[
  rownames(gene_mat_median),
  "gene_name"
]

pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_pec_snRNA_median_logCPM.pdf"
  ),
  height = 20
)
ComplexHeatmap::pheatmap(
  mat = gene_mat_median[
    ,
    order(colnames(gene_mat_median))
  ],
  name = "Scaled median logCPM",
  color = viridis(100, option = "magma"),
  scale = "row",
  column_title = "PEC snRNA",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  row_split = ann_df[rownames(gene_mat_median), ],
  annotation_row = ann_df[rownames(gene_mat_median), , drop = FALSE],
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
