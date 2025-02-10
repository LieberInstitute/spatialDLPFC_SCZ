# Load library ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(here)
  library(spatialLIBD)
  library(ComplexHeatmap)
  library(ggplot2)
  library(sessioninfo)
})

# Load data ----
ntc_enrich_res <- readRDS(here(
  "processed-data/rds/layer_enrich_test",
  "test_enrich_PRECAST_07_ntc.rds"
))

scz_enrich_res <- readRDS(here(
  "processed-data/rds/layer_enrich_test",
  "test_enrich_PRECAST_07_scz.rds"
))


# Make scatter plot ----
t_stat_df <- full_join(
  ntc_enrich_res |> select(ensembl, starts_with("t_stat_")),
  scz_enrich_res |> select(ensembl, starts_with("t_stat_")),
  by = "ensembl",
  suffix = c("_ntc", "_scz")
)


ggplot(t_stat_df) +
  geom_point(aes(x = t_stat_spd01_scz, y = t_stat_spd01_ntc))

with(t_stat_df, cor(t_stat_spd01_scz, t_stat_spd01_ntc))

# Spatial registration plot ----
ntc_t_stats <- ntc_enrich_res[, grep("^t_stat_", colnames(ntc_enrich_res))]
colnames(ntc_t_stats) <- gsub("^t_stat_", "", colnames(ntc_t_stats))


manual_cor_all <- layer_stat_cor(
  ntc_t_stats,
  list("enrichment" = scz_enrich_res),
  model_type = "enrichment",
  reverse = FALSE,
  top_n = 100
)

# Calculate confidence
annotate_registered_clusters(
  cor_stats_layer = manual_cor_all,
  confidence_threshold = 0.45,
  cutoff_merge_ratio = 0.25
)

#   cluster layer_confidence layer_label
# 1   spd01             good spd01/spd04
# 2   spd04             good spd04/spd01
# 3   spd07             good       spd07
# 4   spd02             good       spd02
# 5   spd06             good       spd06
# 6   spd03             good       spd03
# 7   spd05             good       spd05


# Order based on layers
manual_cor_all <- manual_cor_all[sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4)), sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4))]


library(ComplexHeatmap)

my.col <-
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(100)


scz_color_bar <- rowAnnotation(
  ` ` = rownames(manual_cor_all),
  col = list(` ` = set_names(
    Polychrome::palette36.colors(7)[seq.int(7)],
    rownames(manual_cor_all) |> sort()
  )),
  show_legend = FALSE
)

ntc_color_bar <- columnAnnotation(
  ` ` = colnames(manual_cor_all),
  col = list(` ` = set_names(
    Polychrome::palette36.colors(7)[seq.int(7)],
    colnames(manual_cor_all) |> sort()
  )),
  show_legend = FALSE
)

pdf(
  here(
    "plots/05_spatial_registration",
    "heatmap_registration_between_dx_100G.pdf"
  ),
  height = 4,
  width = 5.5
)
Heatmap(
  manual_cor_all,
  name = "Corr",
  col = my.col,
  # row_split = layer_anno_all$bayesSpace,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  right_annotation = scz_color_bar,
  bottom_annotation = ntc_color_bar,
  row_title = "SCZ",
  column_title = "NTC"
  # column_order = ,
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(anno_matrix[i, j], x, y, gp = gpar(fontsize = 10))
  # }
) |> print()
dev.off()

# Session info ----
session_info()