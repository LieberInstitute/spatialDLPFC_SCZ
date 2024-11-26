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
## Load SPG spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

## Load Spatial Domain data ----
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

## Call SPG + spots ----
# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)


# Keep only neun spot =
spe <- spe[ ,spe$neun_pos == TRUE]
# 55,017 spots
stopifnot(ncol(spe) < 80000) # Error prevention


## Create domain-specific neun_claudin label ----
spe$spd_neun <- spe$PRECAST_07






# Create Cell Chat Objects for dx groups ----
## Function ----
subset_N_cellchat <- function(spe, .dx = "ntc", var) {
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
    group.by = var
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

## Domain-agnostic CellChat object ----
set.seed(20241107)
ntc_cellchat <- subset_N_cellchat(spe, .dx = "ntc", var = "spd_neun")
ntc_cellchat |> saveRDS(
  #TODO: change the path
  here(
    "processed-data/rds/spg_ccc",
    "neun_only_ntc_cellchat.rds"
  )
)

rm(ntc_cellchat) # Clean Memory

set.seed(20241107)
scz_cellchat <- subset_N_cellchat(spe, .dx = "scz", var = "spd_neun")
scz_cellchat |> saveRDS(
  here(
    "processed-data/rds/spg_ccc",
    "neun_only_scz_cellchat.rds"
  )
)

# ntc_cellchat <- readRDS(here(
#   "processed-data/rds/spg_ccc",
#   "neun_claudin_ntc_cellchat.rds"
# ))
# scz_cellchat <- readRDS(
#   here(
#     "processed-data/rds/spg_ccc",
#     "neun_claudin_scz_cellchat.rds"
#   )
# )

object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat,
  here(
    "processed-data/rds/spg_ccc",
    "neun_only_merged_cellchat.rds"
  )
)


# Session Info ----
session_info()