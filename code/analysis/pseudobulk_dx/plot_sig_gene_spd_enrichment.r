# Load packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  # library(pheatmap)
  library(ComplexHeatmap)
  library(viridis)
  library(sessioninfo)
})


# Load Data ----
## All data ---
spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

## NTC only ---

ntc_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07_ntc.rds"
  )
)

## SCZ only ---

scz_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

# TODO: remove this
# This matrix isn't that useful actually
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


# Heat map ----

col_df <- colData(spe_pb) |>
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

gene_mat <- logcounts(spe_pb)[sig_gene_df$ensembl, rownames(col_df)]

# Need to accumulate genes
gene_mat_long <- gene_mat |>
  data.frame() |>
  rownames_to_column(var = "gene") |>
  pivot_longer(
    cols = starts_with("V"),
    names_to = c("sample_id", "spd"),
    names_pattern = "(.*)_(spd.*)"
  )

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
  data.matrix()

# Convert Ensembl to Gene name
rownames(gene_mat_median) <- rowData(spe_pb)[
  rownames(gene_mat_median),
  "gene_name"
]


# Format SpD to readable form
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

colnames(gene_mat_median) <- spd_anno_df$anno_lab[
  match(colnames(gene_mat_median), spd_anno_df$spd)
]


gene_mat_median_scaled <- apply(
  gene_mat_median,
  MARGIN = 1,
  FUN = scale
) |> t()

# hc <- hclust(dist(gene_mat_median_scaled))

# gene_names_hc_ordered <- rownames(gene_mat_median)[hc$order]

# saveRDS(
#   gene_names_hc_ordered,
#   here(
#     "code/analysis/pseudobulk_dx",
#     "spd_hierarchical_cluster_order.rds"
#   )
# )


pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_SCZ_spd_median_logCPM.pdf"
  ),
  height = 20
)
heatmap_res <- ComplexHeatmap::pheatmap(
  mat = gene_mat_median[
    # gene_names_hc_ordered,
    ,
    order(colnames(gene_mat_median))
  ],
  name = "Scaled median logCPM",
  color = viridis(100, option = "magma"),
  scale = "row",
  # cluster_rows = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  row_split = ann_df[rownames(gene_mat_median), ],
  annotation_row = ann_df[rownames(gene_mat_median), , drop = FALSE],
  show_row_dend = FALSE,
  # annotation_col = col_df |> select(
  #   PRECAST_07,
  #   # sample_id, # Overwhelm color pallete
  #   dx
  # ),
  cellwidth = 10,
  cellheight = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_colors = list(SCZ_reg = c("Up" = "red", "Down" = "blue"))
)

plot(heatmap_res)

dev.off()


gene_names_hc_ordered <- heatmap_res@row_names_param$labels

saveRDS(
  gene_names_hc_ordered,
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)




# plot(pheatmap_result)

# # Orignial scale ----
# pheatmap(
#   mat = gene_mat,
#   scale = "none",
#   cluster_cols = FALSE,
#   annotation_col = col_df |> select(
#     PRECAST_07,
#     # sample_id, # Overwhelm color pallete
#     dx
#   ),
#   # Turn off names
#   show_rownames = FALSE,
#   show_colnames = FALSE,
# )

# # Row Centered plot -----
# centered_gene_mat <- apply(gene_mat, 1, scale, scale = FALSE) |> t()

# colnames(centered_gene_mat) <- colnames(gene_mat)

# pheatmap(
#   mat = centered_gene_mat,
#   scale = "none",
#   cluster_cols = FALSE,
#   annotation_col = col_df |> select(
#     PRECAST_07,
#     # sample_id, # Overwhelm color pallete
#     dx
#   ),
#   # Turn off names
#   show_rownames = FALSE,
#   show_colnames = FALSE,
# )

# # Row centered and scaled ----
# pheatmap(
#   mat = gene_mat,
#   scale = "row",
#   cluster_cols = FALSE,
#   annotation_col = col_df |> select(
#     PRECAST_07,
#     # sample_id, # Overwhelm color pallete
#     dx
#   ),
#   # Turn off names
#   show_rownames = FALSE,
#   show_colnames = FALSE,
# )

# Session Info ----
session_info()
