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

# Create Cell Chat Objects for dx groups ----
subset_N_cellchat <- function(spe, .sample_id) {
  sub_spe <- spe[, spe$sample_id == .sample_id]
  data.input <- logcounts(sub_spe)
  # Convert emsembl id to gene symbol
  rownames(data.input) <- rowData(sub_spe)$gene_name # Very important step!
  meta <- as.data.frame(SingleCellExperiment::colData(sub_spe))
  meta$samples <- meta$sample_id |> factor()

  # browser()

  # create object
  cellchat <- createCellChat(
    object = data.input, meta = meta,
    group.by = "PRECAST_07"
  )

  # select data base
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  # future::plan("multisession", workers = 3) # do parallel
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
# TODO: create cellchat object for each sample
# Implement a for loop

unique(spe$sample_id) |>
  walk(
    .f = function(.id) {
      cellchat <- subset_N_cellchat(spe, .sample_id = .id)

      saveRDS(
        cellchat,
        here(
          "processed-data/layer_layer_comm/per_sample",
          paste0(.id, ".rds")
        )
      )
    }
  )



# Session Info ----
session_info()
