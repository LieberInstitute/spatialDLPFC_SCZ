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

# Change label ----
# Find annotation file
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") ")) |>
  column_to_rownames("spd")


# Create a column with annoated labels.
spe$spd_annotated <- spd_anno_df[spe$PRECAST_07, "anno_lab"]


# Create Cell Chat Objects for dx groups ----
subset_N_cellchat <- function(spe, .dx = "ntc") {
  sub_spe <- spe[, spe$dx == .dx]
  data.input <- logcounts(sub_spe)
  # Convert emsembl id to gene symbol
  rownames(data.input) <- rowData(sub_spe)$gene_name # Very important step!
  meta <- as.data.frame(SingleCellExperiment::colData(sub_spe))
  meta$samples <- meta$sample_id |> factor()

  # browser()

  # create object
  cellchat <- createCellChat(
    object = data.input, meta = meta,
    group.by = "spd_annotated"
  )

  # select data base
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 3) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  # to prevent memory limit error
  options(future.globals.maxSize = 8.5 * 1024 * 1024^2)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}

set.seed(20241004)
ntc_cellchat <- subset_N_cellchat(spe, .dx = "ntc")
ntc_cellchat |> saveRDS(
  here(
    "processed-data/layer_layer_comm",
    "ntc_cellchat_annotated.rds"
  )
)

rm(ntc_cellchat) # Clean Memory

set.seed(20241004)
scz_cellchat <- subset_N_cellchat(spe, .dx = "scz")
scz_cellchat |> saveRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_annotated.rds"
  )
)


#

# Create merged cellChat object ----
ntc_cellchat <- readRDS(here(
  "processed-data/layer_layer_comm",
  "ntc_cellcha_annotatedt.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_annotated.rds"
  )
)

object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat,
  here(
    "processed-data/layer_layer_comm",
    "merged_cellchat_annotated.rds"
  )
)

# # Compare the total number of interactions & Strength--
# gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
# gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
# gg1 + gg2

# netVisual_diffInteraction(cellchat, weight.scale = T)


# Session Info ----
session_info()
