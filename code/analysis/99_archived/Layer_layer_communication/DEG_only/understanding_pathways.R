library(tidyverse)
library(here)
library(CellChat)
library(circlize)
library(ComplexHeatmap)

# Load data ----
## LR pair enriched in CellChatDB ----
LR_pair_both <- read_csv(
  here(
    "code/analysis/Layer_layer_communication/DEG_only",
    "LR_pair_DEG_both.csv"
  )
)

unique(LR_pair_both$pathway_name)

LR_pair_single <- read_csv(
  here(
    "code/analysis/Layer_layer_communication/DEG_only",
    "LR_pair_DEG_single.csv"
  )
)
unique(LR_pair_single$pathway_name)
unique(LR_pair_single$pathway_name) |> length()

## Load (data driven) layer-layer comm results ----
cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "merged_cellchat_annotated.rds"
  )
)

cellchat_pathways <- rankNet(cellchat, return.data = TRUE)$signaling.contribution


# Pathway-level analysis ----

## Find overlapping Pathways ----
cellchat_pathways |> View()
cellchat_diff_pathways <- unique(cellchat_pathways$name) |> as.character()

cellchat_diff_pathways |> length()

overlap_pathways <- intersect(cellchat_diff_pathways, LR_pair_single$pathway_name)

## Layer-specificity of the pathways ----
# NOTE: the netVisual_heatmap might have bug where it just shows everything
# for (.path in overlap_pathways) {
#   # BUG (CAUTIOUS): Why the output figures are exactly the same
#   # Need more test:
#   netVisual_heatmap(cellchat,
#     signaling = .path,
#     measure = "weight",
#     title.name = paste0(.path, " (Weight)")
#   ) |> print()
# }

ntc_cellchat <- readRDS(here(
  "processed-data/layer_layer_comm",
  "ntc_cellchat_annotated.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_annotated.rds"
  )
)

object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)

# Source code to create diff circular plot
source(here(
  "code/analysis/SPG_analysis/ccc_neun_claudin",
  "3_A-Sender_Receiver_Chord.r"
))
## Heatmap show diff
source(
  here(
    "code/analysis/Layer_layer_communication/vulnerablility_charact",
    "fun_extract_pathway_diff_mat.r"
  )
)

plot_pathway_panel <- function(vec_pathway, plot_path) {
  # TODO: create pdf file to save
  pdf(
    plot_path,
    height = 7, width = 9
  )
  for (.path in vec_pathway) {
    ## 2 circular plots for NTC and SCZ respectively

    par(mfrow = c(1, 2), xpd = TRUE)
    pathways.show <- .path
    weight.max <- getMaxWeight(object.list,
      slot.name = c("netP"),
      attribute = pathways.show
    ) # control the edge weig
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]],
        signaling = pathways.show, layout = "chord",
        edge.weight.max = weight.max[1], edge.width.max = 10,
        signaling.name = paste(pathways.show, names(object.list)[i])
      )
    }
    par(mfrow = c(1, 1), xpd = TRUE)
    ## Circular plots for diff
    my_netVisual_aggregate(ntc_cellchat, scz_cellchat,
      signaling = .path,
      measure = "weight", layout = "chord"
    )

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
      hist = anno_barplot(pathway_matrix |> colSums(), which = "column")
    )


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
    ) |> print()
  }

  dev.off()
}

# DEG overlap pathways ----
plot_pathway_panel(
  vec_pathway = overlap_pathways,
  plot_path = here(
    "plots/layer_layer_comm/",
    "DEG_overlap_pathways_annotated_w_diff.pdf"
  )
)
# DEG overlap pathways ----

# Likeli_pathway Curated by Sang Ho ----
likeli_pathway <- c(
  "NRXN", "L1CAM", "CNTN", "PTN", "GABA-A", "GABA-B",
  "Glutamate"
)

plot_pathway_panel(
  vec_pathway = likeli_pathway,
  plot_path = here(
    "plots/layer_layer_comm/",
    "SCZ_likely_pathways_annotated_w_diff.pdf"
  )
# )

# Potential_pathway Curated by Sang Ho ----
novel_pathway <- c(
  "APP", "CypA", "SPP1")
  # , "NT", "SOMATOSTATIN")
plot_pathway_panel(
  vec_pathway = novel_pathway,
  plot_path =
here(
    "plots/layer_layer_comm/",
    "SCZ_pop_novel_pathways_annotated_w_diff.pdf"
  )
)

# pdf(
#   here(
#     "plots/layer_layer_comm/",
#     "SCZ_likely_pathways_annotated.pdf"
#   ),
#   height = 7, width = 9
# )

# NT and SOMATOSTATIN Pathway
# No NT pathway detected
par(mfrow = c(1, 2), xpd = TRUE)
for (.path in c("NT", "SOMATOSTATIN")) {
  browser()
  pathways.show <- .path
  weight.max <- getMaxWeight(object.list,
    slot.name = c("netP"),
    attribute = pathways.show
  ) # control the edge weig
  for (i in 1:length(object.list)) {
    # browser()
    netVisual_aggregate(object.list[[i]],
      signaling = pathways.show, layout = "chord",
      edge.weight.max = weight.max[1], edge.width.max = 10,
      signaling.name = paste(pathways.show, names(object.list)[i])
    )
  }
}
par(mfrow = c(1, 1), xpd = TRUE)
# dev.off()

# TODO: create pdf file to save
# pdf(
#   here(
#     "plots/layer_layer_comm/",
#     "SCZ_pop_novel_pathways_annotated.pdf"
#   ),
#   height = 7, width = 9
# )

# par(mfrow = c(1, 2), xpd = TRUE)
# for (.path in novel_pathway) {
#   pathways.show <- .path
#   weight.max <- getMaxWeight(object.list,
#     slot.name = c("netP"),
#     attribute = pathways.show
#   ) # control the edge weig
#   for (i in 1:length(object.list)) {
#     # browser()
#     netVisual_aggregate(object.list[[i]],
#       signaling = pathways.show, layout = "chord",
#       edge.weight.max = weight.max[1], edge.width.max = 10,
#       signaling.name = paste(pathways.show, names(object.list)[i])
#     )
#   }
# }
# par(mfrow = c(1, 1), xpd = TRUE)
# dev.off()

# TODO: possibly plotting the difference plot


# netVisual_diffInteraction

# pdf(
#   here(
#     "plots/layer_layer_comm/",
#     "SCZ_likely_pathways_diff.pdf"
#   )
# )
# for (.path in likeli_pathway) {
#   print(my_netVisual_aggregate(ntc_cellchat, scz_cellchat, signaling = .path, measure = "weight", layout = "chord"))
# }
# dev.off()




## Where pathway differs ----




# LR level analysis ----
# How many LR pairs ----

# overlapping pairs  ---
# between enrichment and data-driven layer-layer analysis

# What are those

# Layer specificity

# Differential
