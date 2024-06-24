# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(pheatmap)
  library(sessioninfo)
})

# Read Data ----
## Dx DEG genes ----
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

ann_df <- sig_gene_df  |>
  # filter(fdr_ntc <= 0.05) |>
  column_to_rownames(var = "gene") |>
  transmute(
    SCZ_reg = factor(
      logFC_scz > 0,
      levels = c(TRUE, FALSE),
      labels = c("Up", "Down")
    )
  )

spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

# Enrichment test results (-log10(p-value)) ----
modeling_results <- readRDS(
  here(
    "processed-data", "rds", "layer_enrich_test",
    "test_enrich_PRECAST_07.rds"
  )
)

rownames(modeling_results) <- NULL

heatmap_pec_spd_df <- modeling_results |>
  mutate(
    across(
      starts_with("p_value"),
      ~ -1 * log10(.x),
      .names = "neg_log10_{.col}"
    )
  ) |>
  filter(ensembl %in% sig_gene_df$ensembl) |>
  column_to_rownames(var = "gene") |>
  select(starts_with("neg_log10_")) |>
  rename_with(
    .fn = ~ gsub("^neg_log10_p_value_", "", .),
    .cols = starts_with("neg_log10_")
  )

colnames(heatmap_pec_spd_df) <- spd_anno_df$anno_lab[
  match(colnames(heatmap_pec_spd_df), spd_anno_df$spd)
]


hc <- hclust(dist(heatmap_pec_spd_df))

gene_names_hc_ordered <- rownames(heatmap_pec_spd_df)[hc$order]

saveRDS(
  gene_names_hc_ordered,
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)



pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_SCZ_spd_p_value.pdf"
  ),
  height = 20
)
heatmap_pec_spd_df[gene_names_hc_ordered, order(colnames(heatmap_pec_spd_df))] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    # cluster_rows = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )
dev.off()

# Plot t-testistics ----

t_stat_mat <- modeling_results |>
  mutate(
    across(
      starts_with("t_stat"),
      ~ abs(.x),
      .names = "abs_{.col}"
    )
  ) |>
  filter(ensembl %in% sig_gene_df$ensembl) |>
  column_to_rownames(var = "gene") |>
  select(starts_with("abs_t_stat_")) |>
  rename_with(
    .fn = ~ gsub("^abs_t_stat_", "", .),
    .cols = starts_with("abs_t_")
  )

colnames(t_stat_mat) <- spd_anno_df$anno_lab[
  match(colnames(t_stat_mat), spd_anno_df$spd)
]


pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_SCZ_spd_t_stat.pdf"
  ),
  height = 20
)
t_stat_mat[hc$order, order(colnames(t_stat_mat))] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    # scale = "row",
    # cluster_rows = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )
dev.off()

pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_SCZ_spd_t_stat_row_centered.pdf"
  ),
  height = 20
)
t_stat_mat[hc$order, order(colnames(t_stat_mat))] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    # cluster_rows = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )

dev.off()




# Session Info ----
session_info()
