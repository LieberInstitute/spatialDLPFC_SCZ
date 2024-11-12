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

## Create variable for neun_claudin ----
# TODO: this ineed to be edtited to define those both spots
spe$neun_claudin <- case_when(
  # NOTE: if a spot is both neun and claudin +, it is defnied as neun+ spot
  spe$neun_pos == TRUE  ~ "neun",
  spe$vasc_pos == TRUE ~ "claudin",
  .default = NA_character_
)

# Remove spots that are not in either of these category
spe <- spe[,!is.na(spe$neun_claudin)]
stopifnot(ncol(spe) < 80000) # Error prevention



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
    group.by = "neun_claudin"
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

set.seed(20241107)
ntc_cellchat <- subset_N_cellchat(spe, .dx = "ntc")
ntc_cellchat |> saveRDS(
  #TODO: change the path
  here(
    "processed-data/rds/spg_ccc",
    "neun_claudin_ntc_cellchat.rds"
  )
)

rm(ntc_cellchat) # Clean Memory

set.seed(20241107)
scz_cellchat <- subset_N_cellchat(spe, .dx = "scz")
scz_cellchat |> saveRDS(
  here(
    "processed-data/rds/spg_ccc",
    "neun_claudin_scz_cellchat.rds"
  )
)


#

# Create merged cellChat object ----
ntc_cellchat <- readRDS(here(
  "processed-data/rds/spg_ccc",
  "neun_claudin_ntc_cellchat.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/rds/spg_ccc",
    "neun_claudin_scz_cellchat.rds"
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
    "processed-data/rds/spg_ccc",
    "neun_claudin_merged_cellchat.rds"
  )
)

# Session Info ----
session_info()