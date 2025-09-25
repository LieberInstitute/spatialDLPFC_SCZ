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

# Subset DEGs only ----
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
)

sig_genes <- gene_df |> filter(fdr_scz <= 0.10)
# sig_genes$ensembl

# spe <- spe[sig_genes$ensembl, ]

# stopifnot(nrow(spe) == nrow(sig_genes))


# Create Cell Chat Objects for dx groups ----
subset_N_cellchat <- function(spe, .dx = "ntc") {
  sub_spe <- spe[, spe$dx == .dx]
  data.input <- logcounts(sub_spe)
  # Convert emsembl id to gene symbol
  rownames(data.input) <- rowData(sub_spe)$gene_name # Very important step!
  meta <- as.data.frame(SingleCellExperiment::colData(sub_spe))
  meta$samples <- meta$sample_id |> factor()



  # create object
  cellchat <- createCellChat(
    object = data.input, meta = meta,
    group.by = "PRECAST_07"
  )
  browser()
  # select data base
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  # future::plan("multisession", workers = 3) # do parallel
  # I don't feel like we need to run this step, but don't know why we have to.
  cellchat <- identifyOverExpressedGenes(cellchat)
  # cellchat@meta

  # Using DEGs in the CellChat pipeline.
  cellchat@var.features[["features"]] <- sig_genes$gene
  # to prevent memory limit error
  options(future.globals.maxSize = 8.5 * 1024 * 1024^2)

  write_csv(
    identifyOverExpressedInteractions(
      cellchat,
      features = sig_genes$gene,
      return.object = FALSE, variable.both = TRUE
    ),
    here(
      "code/analysis/Layer_layer_communication/DEG_only",
      "LR_pair_DEG_both.csv"
    )
  )

  write_csv(
    identifyOverExpressedInteractions(
      cellchat,
      features = sig_genes$gene,
      return.object = FALSE, variable.both = FALSE
    ),
    here(
      "code/analysis/Layer_layer_communication/DEG_only",
      "LR_pair_DEG_single.csv"
    )
  )

  cellchat <- identifyOverExpressedInteractions(
    cellchat,
    features = sig_genes$gene,
    return.object = TRUE, variable.both = FALSE
  )

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
    "ntc_cellchat_deg_only.rds"
  )
)

rm(ntc_cellchat) # Clean Memory

set.seed(20241004)
scz_cellchat <- subset_N_cellchat(spe, .dx = "scz")
scz_cellchat |> saveRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_deg_only.rds"
  )
)


#

# Create merged cellChat object ----
ntc_cellchat <- readRDS(here(
  "processed-data/layer_layer_comm",
  "ntc_cellchat_deg_only.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_deg_only.rds"
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
    "merged_cellchat_deg_only.rds"
  )
)

# Compare the total number of interactions & Strength--
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
gg1 + gg2

netVisual_diffInteraction(cellchat, weight.scale = T)


# Session Info ----
session_info()
