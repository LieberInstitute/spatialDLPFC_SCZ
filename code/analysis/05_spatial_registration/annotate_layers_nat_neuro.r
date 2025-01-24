# Load library ----
suppressPackageStartupMessages({
  library(spatialLIBD)
  library(here)
  library(ComplexHeatmap)
  library(tidyverse)
  library(sessioninfo)
})


# Load data -----
## Enrichment analysis -----
layer_res <- readRDS(here(
  "processed-data", "rds", "layer_enrich_test",
  paste0("test_enrich_PRECAST_07.rds")
))

## format enrichment test res
t_stats <- layer_res[, grep("^t_stat_", colnames(layer_res))]
colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))

## nat neuro enrichment res ----
# note: doesn't work due to data fetching reason
# manual_modeling_results <- fetch_data(type = "modeling_results")
load(
  here(
    "code/analysis/05_spatial_registration",
    "Human_DLPFC_Visium_modeling_results.Rdata"
  )
)
manual_modeling_results <- modeling_results

# Calculate correlation mat ----
manual_cor_all <- layer_stat_cor(
  t_stats,
  manual_modeling_results,
  model_type = "enrichment",
  reverse = FALSE,
  top_n = 100
)

# Annotation confidence -----
anno_conf_df <- annotate_registered_clusters(
  cor_stats_layer = manual_cor_all,
  confidence_threshold = 0.25,
  cutoff_merge_ratio = 0.25
) |>
  arrange(cluster)
# r$> annotate_registered_clusters(
#       cor_stats_layer = manual_cor_all,
#       confidence_threshold = 0.25,
#       cutoff_merge_ratio = 0.1
#     ) |> arrange(cluster)
#   cluster layer_confidence layer_label
# 1   spd01             good          WM
# 2   spd02             good        L4/3
# 3   spd03             good          L6
# 4   spd04             good          WM
# 5   spd05             good          L5
# 6   spd06             good        L2/3
# 7   spd07             good          L1

anno_mat_long <- anno_conf_df |>
  pmap_dfr(.f = function(cluster, layer_label, ...) {
    # browser()
    data.frame(
      cluster = cluster,
      layer_label = str_split(layer_label, "/") |> unlist()
    )
  }) |>
  mutate(
    layer_label = case_when(
      str_length(layer_label) == 1 ~ paste0("L", layer_label),
      TRUE ~ layer_label
    )
  ) |>
  mutate(
    layer_label = str_replace(layer_label, "L", "Layer")
  )


# Heatmap ----
my.col <-
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(100)

layer_color_bar <- columnAnnotation(
  `Manual\nAnnotation` = colnames(manual_cor_all),
  col = list(`Manual\nAnnotation` = spatialLIBD::libd_layer_colors),
  show_legend = FALSE,
  annotation_name_side = "left"
)

subset_color_bar <- rowAnnotation(
  ` ` = rownames(manual_cor_all),
  col = list(` ` = set_names(
    Polychrome::palette36.colors(7)[seq.int(7)],
    rownames(manual_cor_all) |> sort()
  )),
  show_legend = FALSE
)


anno_matrix <- matrix(" ",
  nrow = nrow(manual_cor_all),
  ncol = ncol(manual_cor_all)
)

rownames(anno_matrix) <- rownames(manual_cor_all)
colnames(anno_matrix) <- colnames(manual_cor_all)


n_good <- nrow(anno_mat_long)
for (i in 1:n_good) {
  anno_matrix[
    anno_mat_long$cluster[i],
    anno_mat_long$layer_label[i]
  ] <- "X"
}

anno_matrix <- anno_matrix[
  anno_conf_df$cluster[order(anno_conf_df$layer_label)],
  sort(colnames(manual_cor_all))
]

pdf(
  here("plots/05_spatial_registration",
  # "heatmap_registration_nat_neuro_100G_with_X.pdf"),
    "heatmap_registration_nat_neuro_100G.pdf"),
  height = 4,
  width = 5.5
)
Heatmap(
  manual_cor_all[
    anno_conf_df$cluster[order(anno_conf_df$layer_label)],
    sort(colnames(manual_cor_all))
  ],
  name = "Cor",
  col = my.col,
  # row_split = layer_anno_all$bayesSpace,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  right_annotation = subset_color_bar,
  bottom_annotation = layer_color_bar,
  # column_order = ,
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(anno_matrix[i, j], x, y, gp = gpar(fontsize = 10))
  # }
) |> print()
dev.off()

# Save 



# session Info ----
session_info()
