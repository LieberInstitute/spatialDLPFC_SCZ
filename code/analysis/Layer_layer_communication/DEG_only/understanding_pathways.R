library(tidyverse)
library(here)
library(CellChat)

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
    "merged_cellchat.rds"
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
  "ntc_cellchat.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat.rds"
  )
)

object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)

# TODO: create pdf file to save
pdf(
  here(
    "plots/layer_layer_comm/",
    "DEG_overlap_pathways_annotated.pdf"
  ),
  height = 7, width = 9
)
par(mfrow = c(1, 2), xpd = TRUE)
for (.path in overlap_pathways) {
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
dev.off()


# TODO: create pdf file to save
likeli_pathway <- c(
  "NRXN", "L1CAM", "CNTN", "PTN", "GABA-A", "GABA-B",
  "Glutamate"
)

pdf(
  here(
    "plots/layer_layer_comm/",
    "SCZ_likely_pathways_annotated.pdf"
  ),
  height = 7, width = 9
)
par(mfrow = c(1, 2), xpd = TRUE)
for (.path in likeli_pathway) {
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
dev.off()

# TODO: create pdf file to save
pdf(
  here(
    "plots/layer_layer_comm/",
    "SCZ_pop_novel_pathways_annotated.pdf"
  ),
  height = 7, width = 9
)
novel_pathway <- c(
  "CypA", "SPP1", "APP"
)
par(mfrow = c(1, 2), xpd = TRUE)
for (.path in novel_pathway) {
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
dev.off()

# TODO: possibly plotting the difference plot


# netVisual_diffInteraction
source(here(
  "code/analysis/SPG_analysis/ccc_neun_claudin",
  "3_A-Sender_Receiver_Chord.r"
))

pdf(
  here(
    "plots/layer_layer_comm/",
    "SCZ_likely_pathways_diff.pdf"
  )
)
for (.path in likeli_pathway) {
  print(my_netVisual_aggregate(ntc_cellchat, scz_cellchat, signaling = .path, measure = "weight", layout = "chord"))
}
dev.off()




## Where pathway differs ----




# LR level analysis ----
# How many LR pairs ----

# overlapping pairs  ---
# between enrichment and data-driven layer-layer analysis

# What are those

# Layer specificity

# Differential
