# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(sessioninfo)
})

# Load data ----
raw_chea_df <-
  read_csv(
    here("processed-data/rds/13_TF_enrichment", "ChEA3_result_compiled.csv")
  )

up_chea_df <- raw_chea_df |>
  filter(direction == "up")

down_chea_df <- raw_chea_df |>
  filter(direction == "down")


up_chea_df_wide <- up_chea_df |>
  pivot_wider(id_cols = TF, names_from = spd, values_from = Rank) |>
  select(c(
    "TF", "SpD07_L1", "SpD06_L2_3", "SpD02_L3_4",
    "SpD05_L5", "SpD03_L6", "SpD01_WMtz", "SpD04_WM"
  ))

down_chea_df_wide <- down_chea_df |>
  pivot_wider(id_cols = TF, names_from = spd, values_from = Rank) |>
  select(c(
    "TF", "SpD07_L1", "SpD06_L2_3", "SpD02_L3_4",
    "SpD05_L5", "SpD03_L6", "SpD01_WMtz", "SpD04_WM"
  ))

n_tf_total <- nrow(up_chea_df_wide)
stopifnot(
  n_tf_total == nrow(down_chea_df_wide)
)

spd_name_df <-
  data.frame(
    raw = c(
      "SpD07_L1", "SpD06_L2_3", "SpD02_L3_4",
      "SpD05_L5", "SpD03_L6", "SpD01_WMtz", "SpD04_WM"
    ),
    final = c(
      "SpD07-L1/M", "SpD06-L2/3", "SpD02-L3/4",
      "SpD05-L5", "SpD03-L6", "SpD01-WMtz", "SpD04-WM"
    )
  )

# Identify overlap -----
## GWAS genes (trubetskoy et al) ----
# NOTE: supp table 12 from Trubetskoy
trubetskoy_df <- read_excel(
  here(
    "code/analysis/Layer_layer_communication/genetic_risk",
    "Supplementary Table 12.xlsx"
  ),
  sheet = "Prioritised"
)

symbols_gene <- trubetskoy_df |>
  # filter(gene_biotype == "protein_coding") |>
  pull(Symbol.ID)

# Trubetskoy Overlap with upregulated TFs ----
up_chea_df_wide |>
  filter(TF %in% c(symbols_gene)) |>
  View()



# TF enriched in upreg DEGs ----
org_up_trub_tf_wide <- up_chea_df_wide |>
  filter(TF %in% c(symbols_gene)) |>
  column_to_rownames("TF") |>
  as.matrix()
org_up_trub_tf_wide <- org_up_trub_tf_wide / n_tf_total

colnames(org_up_trub_tf_wide) <- spd_name_df$final[match(colnames(org_up_trub_tf_wide), spd_name_df$raw)]


p_trub_tf_up <- Heatmap(
  matrix = org_up_trub_tf_wide |> t(), # Normalize by total number of TFs
  name = "Up-reg ",
  col = colorRamp2(
    c(0, 0.5, 1), # Reverse the scale: low rank = black, high rank = white
    c("red", "white", "blue")
  ),
  cluster_rows = FALSE, # Same SpD orders
  cluster_columns = TRUE, # Cluster the TFs
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_title = NULL, # Do not show row data titles
  column_title = NULL, # Do not show row
  heatmap_legend_param = list(
    title = "Rank Percentile",
    at = c(0, 0.5, 1),
    labels = c(
      "highly ranked (relavent)",
      "Median",
      "lowly ranked (not relevant)"
    ),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8)
  ),
  row_names_gp = gpar(fontsize = 8, col = "red"),
  column_names_gp = gpar(fontsize = 8),
  show_row_dend = FALSE,
  show_column_dend = FALSE # Hide  # Hide dendrogram but still cluster rows
)

pdf(
  here("plots/13_TF_enrichment/Trubetskoy_TF_enrichment_upregulated_DEGs.pdf"),
  height = 1.5, width = 2
)
# Add a title using draw()
draw(p_trub_tf_up,
  # column_title = "Enrichment of Upregulated layer-restricted DEGs",
  # column_title_gp = gpar(fontsize = 8)
  show_heatmap_legend = FALSE
)
dev.off()


# Trubetskoy Overlap with downregulated TFs ----
down_chea_df_wide |>
  filter(TF %in% c(symbols_gene))

org_down_trub_tf_wide <- down_chea_df_wide |>
  filter(TF %in% c(symbols_gene)) |>
  column_to_rownames("TF") |>
  as.matrix()
org_down_trub_tf_wide <- org_down_trub_tf_wide / n_tf_total

colnames(org_down_trub_tf_wide) <- spd_name_df$final[match(colnames(org_down_trub_tf_wide), spd_name_df$raw)]


p_trub_tf_down <- Heatmap(
  matrix = org_down_trub_tf_wide |> t(), # Normalize by total number of TFs
  name = "Down-reg ",
  col = colorRamp2(
    c(0, 0.5, 1), # Reverse the scale: low rank = black, high rank = white
    c("red", "white", "blue")
  ),
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = NULL,
  row_title = NULL,
  heatmap_legend_param = list(
    title = "Rank Percentile",
    at = c(0, 0.5, 1),
    labels = c(
      "highly ranked (relavent)",
      "Median",
      "lowly ranked (not relevant)"
    )
  ),
  row_names_gp = gpar(fontsize = 8, col = "blue"),
  column_names_gp = gpar(fontsize = 8),
  show_row_dend = FALSE,
  show_column_dend = FALSE # Hide  # Hide dendrogram but still cluster rows
)

pdf(
  here("plots/13_TF_enrichment/Trubetskoy_TF_enrichment_downregulated_DEGs.pdf"),
  height = 1.5, width = 2
)

draw(
  p_trub_tf_down,
  show_heatmap_legend = FALSE
  # column_title = "Enrichment of Downregulated layer-restricted DEGs",
  # column_title_gp = gpar(fontsize = 8)
)
dev.off()

# Top ranking TF per domains -----

## TF from up-reg DEGs -----
subset_top_up_tf <- up_chea_df |>
  group_by(spd) |>
  slice_min(Rank, n = 10)

# Number of unique TF among top TFs
subset_top_up_tf |>
  pull(TF) |>
  unique() |>
  length()

org_up_top_tf_wide <- up_chea_df_wide |>
  filter(TF %in%
    (subset_top_up_tf |> pull(TF) |> unique())) |>
  column_to_rownames("TF") |>
  as.matrix()
org_up_top_tf_wide <- org_up_top_tf_wide / n_tf_total

colnames(org_up_top_tf_wide) <- spd_name_df$final[match(colnames(org_up_top_tf_wide), spd_name_df$raw)]

p_top_up_tf <- Heatmap(
  matrix = org_up_top_tf_wide |> t(), # Normalize by total number of TFs
  name = "Up-reg ",
  col = colorRamp2(
    c(0, 0.5, 1), # Reverse the scale: low rank = black, high rank = white
    c("red", "white", "blue")
  ),
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = NULL,
  row_title = NULL,
  heatmap_legend_param = list(
    title = "Rank Percentile",
    at = c(0, 0.5, 1),
    labels = c(
      "highly ranked (relavent)",
      "Median",
      "lowly ranked (not relevant)"
    )
  ),
  row_names_gp = gpar(fontsize = 8, col = "red"),
  column_names_gp = gpar(fontsize = 8),
  show_row_dend = FALSE,
  show_column_dend = FALSE # Hide  # Hide dendrogram but still cluster rows
)

pdf(
  here("plots/13_TF_enrichment/Top_TF_enrichment_upregulated_DEGs.pdf"),
  height = 1.5, width = 4.5
)
draw(p_top_up_tf,
  # column_title = "Enrichment of Top Upregulated layer-restricted DEGs",
  # column_title_gp = gpar(fontsize = 8)
  show_heatmap_legend = FALSE
)
dev.off()

## TF from down-reg DEGs -----
subset_down_tf <- down_chea_df |>
  group_by(spd) |>
  slice_min(Rank, n = 10)

# Number of unique TF among top TFs
subset_down_tf |>
  pull(TF) |>
  unique() |>
  length()

org_down_top_tf_wide <- down_chea_df_wide |>
  filter(TF %in%
    (subset_down_tf |> pull(TF) |> unique())) |>
  column_to_rownames("TF") |>
  as.matrix()
org_down_top_tf_wide <- org_down_top_tf_wide / n_tf_total

colnames(org_down_top_tf_wide) <- spd_name_df$final[match(colnames(org_down_top_tf_wide), spd_name_df$raw)]

p_top_down_tf <- Heatmap(
  matrix = org_down_top_tf_wide |> t(), # Normalize by total number of TFs
  name = "Down-reg ",
  col = colorRamp2(
    c(0, 0.5, 1), # Reverse the scale: low rank = black, high rank = white
    c("red", "white", "blue")
  ),
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = NULL,
  row_title = NULL,
  heatmap_legend_param = list(
    title = "Rank Percentile",
    at = c(0, 0.5, 1),
    labels = c(
      "highly ranked (relavent)",
      "Median",
      "lowly ranked (not relevant)"
    )
  ),
  row_names_gp = gpar(fontsize = 8, col = "blue"),
  column_names_gp = gpar(fontsize = 8),
  show_row_dend = FALSE,
  show_column_dend = FALSE # Hide  # Hide dendrogram but still cluster rows
)

pdf(
  here("plots/13_TF_enrichment/Top_TF_enrichment_downregulated_DEGs.pdf"),
    height = 1.5, width = 4.5
)
draw(p_top_down_tf,
  # column_title = "Enrichment of Top Downregulated layer-restricted DEGs",
  # column_title_gp = gpar(fontsize = 8)
  show_heatmap_legend = FALSE
)
dev.off()


# Session info ----
session_info()
