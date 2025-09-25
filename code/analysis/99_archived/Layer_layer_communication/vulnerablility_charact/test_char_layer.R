# Load libraries ----
library(CellChat)
library(tidyverse)
library(here)
library(sessioninfo)


# Load Data ----
## Load (data driven) layer-layer comm results ----
# cellchat <- readRDS(
#   here(
#     "processed-data/layer_layer_comm",
#     "merged_cellchat_annotated.rds"
#   )
# )

ntc_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "ntc_cellchat_annotated.rds"
  )
)

scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_annotated.rds"
  )
)

# TODO: source the matrix extraction function
source(
  here(
    "code/analysis/Layer_layer_communication/vulnerablility_charact",
    "fun_extract_pathway_diff_mat.r"
  )
)

.path <- "JAM"


# Visualize the pathway matrix using a heatmap
library(pheatmap)
extract_diff_pathway(
  ntc_cellchat, scz_cellchat,
  signaling = .path, measure = "weight"
) |>
as.matrix () |>
library(ComplexHeatmap)

# Extract the pathway matrix
pathway_matrix <- extract_diff_pathway(
  ntc_cellchat, scz_cellchat,
  signaling = .path, measure = "weight"
) |>
as.matrix()

# Create row and column annotations
row_annotation <- rowAnnotation(
  hist = anno_barplot(pathway_matrix |> rowSums(), which = "row")
)
col_annotation <- HeatmapAnnotation(
  hist = anno_barplot(pathway_matrix|> colSums(), which = "column")
)

library(circlize)
# Visualize the pathway matrix using a heatmap with annotations
Heatmap(
  pathway_matrix,
  name = .path,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Sender",
  column_title = "Receiver",
  col = colorRamp2(c(min(pathway_matrix), 0, max(pathway_matrix)), c("blue", "white", "red")),
  top_annotation = col_annotation,
  left_annotation = row_annotation
)


# tmp_obj <- rankNet(cellchat, return.data = TRUE)$signaling.contribution


# Session ----
session_info()
