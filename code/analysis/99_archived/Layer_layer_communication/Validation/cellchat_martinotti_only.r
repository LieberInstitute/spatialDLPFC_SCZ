# srun --pty --x11 --mem=150G --time=5:00:00 --cpus-per-task=4 bash
# Load packages ----
suppressPackageStartupMessages({
  if (!require("CellChat")) {
    devtools::install_github("jinworks/CellChat")
  }
  library(CellChat)
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
  library(patchwork)
  library(sessioninfo)
  if (!require("presto")) {
    devtools::install_github("immunogenomics/presto")
  }
  library(presto)
})



# Load Data ----
## Load Spe Object that has spatial domain label ----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

PRECAST_df_final <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

col_data_df <- PRECAST_df_final |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)


## Subset to NTC sample only ----
spe <- spe[, spe$dx == "ntc"]

## Subset to L5 and L1 sample only -----
# spd07 -> Layer 1, spd05 -> Layer 5
spe <- spe[, spe$PRECAST_07 %in% c("spd07", "spd05")]

## Define L5 Martinotti+ and Martinotti- and L1 -----
spe$martinotti <- {
  logcounts(spe)[which(rowData(spe)$gene_name == "PANTR1"), ] > 0 |
    logcounts(spe)[which(rowData(spe)$gene_name == "RELN"), ] > 0
}

spe$new_cat <- case_when(
  spe$martinotti == TRUE ~ "mart+ (L5)",
  spe$martinotti == FALSE ~ NA_character_,
  TRUE ~ NA_character_
)
spe$new_cat[spe$PRECAST_07 == "spd07"] <- "Layer 1"

# Remove Martinotti - spots ----
spe <- spe[, !is.na(spe$new_cat)]


# Cell Chat ----
## Prepare input data ----
spe$samples <- spe$sample_id |> factor()
data.input <- logcounts(spe)
# Convert emsembl id to gene symbol
rownames(data.input) <- rowData(spe)$gene_name # Very important step!
meta <- as.data.frame(SingleCellExperiment::colData(spe))
cellchat <- createCellChat(
  object = data.input, meta = meta,
  group.by = "new_cat"
)

cellchat@DB <- CellChatDB.human

## Run cell chat ----
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
set.seed(20241018)
# NOTE: I didn't use future as suggested by the tutorial
# That's because the cellchat author doesn't seem to control random seed properly in impletation
# options(future.globals.maxSize = 8.5 * 1024 * 1024^2)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)

saveRDS(
  cellchat,
  here(
    "processed-data/layer_layer_comm",
    "ntc_L1_5_cellchat_martinotti_only.rds"
  )
)

# netVisual_aggregate(cellchat, signaling = "GABA-A")
# netAnalysis_contribution(cellchat, signaling = "GABA-A")
# netVisual_bubble(cellchat, signaling = "GABA-A", remove.isolate = FALSE)
