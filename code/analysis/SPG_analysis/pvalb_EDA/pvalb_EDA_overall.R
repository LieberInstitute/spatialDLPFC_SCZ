# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
  library(ggbeeswarm)
  library(escheR)
  library(sessioninfo)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

## Load SpD data ----
finalized_spd <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Attach SpD label to spe ----
col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

## Call PVALB+ Spots----
pvalb_index <- which(rowData(spe)$gene_name == "PVALB")
spe$pvalb_pos <- logcounts(spe)[pvalb_index, ] > 0
spe$logcount_pvalb <- logcounts(spe)[pvalb_index, ]

# Descriptive analysis ----

## Proportion of Pvalb+ spots ----
### Overall ----
col_df <- colData(spe)
col_df |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(
    n = n(),
    prop_pvalb = sum(pvalb_pos) / n(),
    dx = unique(dx)
  ) |>
  ungroup() |>
  ggplot(aes(x = dx, y = prop_pvalb)) +
  geom_boxplot(aes(group = dx)) +
  geom_jitter(aes(color = dx), alpha = 0.7) +
  theme_minimal()

### Per spatial domain ----
col_df |>
  data.frame() |>
  group_by(sample_id, PRECAST_07) |>
  summarize(
    n = n(),
    prop_pvalb = sum(pvalb_pos) / n(),
    dx = unique(dx)
  ) |>
  ungroup() |>
  ggplot(aes(x = PRECAST_07, y = prop_pvalb, fill = dx)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_minimal() +
  scale_x_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  )


## Gene Expression -----
# Violin plot for the expression of Pvalb for all samples
# Not recommend to run due to number of data points
# scater::plotExpression(
#   spe, rownames(spe)[pvalb_index],
#   x = "sample_id", colour_by = "dx"
# )

col_df |>
  data.frame() |>
  filter(pvalb_pos) |>
  ggplot() +
  geom_boxplot(aes(x = sample_id, y = logcount_pvalb, color = dx)) +
  theme_minimal() +
  scale_x_discrete(
    limits = metadata(spe)$dx$sample_id[order(metadata(spe)$dx$dx)]
  )

## Heatmap showing per spatial domain ----
# NOTE: Here we visualize the median expression of PVALB+ spots per SpD and sample
# This is different from visaulize the mean/median PVALB expression of per SpD-sample_id
col_df |>
  data.frame() |>
  filter(pvalb_pos) |>
  group_by(sample_id, PRECAST_07) |>
  summarize(
    median_pvalb = median(logcount_pvalb),
    dx = unique(dx)
  )
library(ComplexHeatmap)

heatmap_data <- col_df |>
  data.frame() |>
  filter(pvalb_pos) |>
  group_by(sample_id, PRECAST_07) |>
  summarize(
    median_pvalb = median(logcount_pvalb),
    dx = unique(dx)
  ) |>
  pivot_wider(names_from = PRECAST_07, values_from = median_pvalb) |>
  column_to_rownames("sample_id")

# Create row annotation for the heatmap
row_annotation <- rowAnnotation(
  dx = col_df$dx[match(rownames(heatmap_data), col_df$sample_id)],
  col = list(
    dx = c("ntc" = "blue", "scz" = "red")
  ),
  annotation_legend_param = list(dx = list(title = "Diagnosis"))
)

# Generate the heatmap with row annotation
# Doesn't seem to have specific pattern.
Heatmap(
  as.matrix(heatmap_data |> select(-dx)),
  name = "Median PVALB Expression",
  row_title = "Sample ID",
  column_title = "Spatial Domain",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  heatmap_legend_param = list(title = "Expression"),
  left_annotation = row_annotation
)


# Session Info ----
sessioninfo::session_info()
