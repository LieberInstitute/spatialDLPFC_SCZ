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
    "code/analysis/SPG_analysis/ccc_ex_neun/",
    "ex_neun_spe.rds"
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

# Call PNN & PVALB + ----
pvalb_index <- which(rowData(spe)$gene_name == "PVALB")
spe$pvalb_pos <- logcounts(spe)[pvalb_index, ] > 0
spe$logcount_pvalb <- logcounts(spe)[pvalb_index, ]

spe$pnn_pvalb <- spe$pnn_pos & spe$pvalb_pos

# Descriptive statistics ---
# Data only contains neun+ spots
table(spe$neun_pos)
spe

# Extremely limited number of spots are labeled as pnn_pos spots
# like 4k_out of
table(spe$pnn_pos)

# Among the neun+ spots, half of them are classified as excitatory neunrons
table(spe$ex_cell)

table(spe$pnn_pvalb, spe$pnn_pos)

# Check
table(spe$pnn_pos, spe$ex_cell)
table(spe$pnn_pvalb, spe$ex_cell)

# Extremely small number of PNN PVALB spots

# Any pnn pvalb spots is pvalb spots to boost sample size

spe$ex_ccc_cat <- case_when(
  spe$pnn_pvalb == TRUE  ~ "PVALB",
  spe$ex_cell == TRUE ~ "Ex Only",
  # spe$pnn_pvalb == TRUE & spe$ex_cell == FALSE ~ "PVALB Only",
  TRUE ~ NA_character_
)



# spe$ex_ccc_cat <- case_when(
#   spe$pnn_pvalb == TRUE & spe$ex_cell == TRUE ~ "PVALB & Ex",
#   spe$pnn_pvalb == FALSE & spe$ex_cell == TRUE ~ "Ex Only",
#   spe$pnn_pvalb == TRUE & spe$ex_cell == FALSE ~ "PVALB Only",
#   TRUE ~ NA_character_
# )


# Remove spots that are not in either of these category
spe <- spe[, !is.na(spe$ex_ccc_cat)]
# GM region - L1 only
spe <- spe[, !spe$PRECAST_07 %in% paste0("spd", c("01", "04", "07"))]
# stopifnot(ncol(spe) < 80000) # Error prevention



## Create domain-specific neun_claudin label ----
spe$spd_pvalb_ex <- paste0(spe$ex_ccc_cat, "_", spe$PRECAST_07)
spe$spd_pvalb_ex |> table()





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
set.seed(20241202)
ntc_cellchat <- subset_N_cellchat(spe, .dx = "ntc", var = "spd_pvalb_ex")
ntc_cellchat |> saveRDS(
  # TODO: change the path
  here(
    "processed-data/rds/spg_ccc",
    "neun_ex_ntc_cellchat.rds"
  )
)

rm(ntc_cellchat) # Clean Memory

set.seed(20241202)
scz_cellchat <- subset_N_cellchat(spe, .dx = "scz", var = "spd_pvalb_ex")
scz_cellchat |> saveRDS(
  here(
    "processed-data/rds/spg_ccc",
    "neun_ex_scz_cellchat.rds"
  )
)

ntc_cellchat <- readRDS(here(
  "processed-data/rds/spg_ccc",
  "spd_neun_ex_ntc_cellchat.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/rds/spg_ccc",
    "spd_neun_ex_scz_cellchat.rds"
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
    "spd_neun_ex_merged_cellchat.rds"
  )
)

# ## Domain-specific CellChat object ----
# set.seed(20241107)
# ntc_cellchat <- subset_N_cellchat(spe, .dx = "ntc", var = "spd_neun_claudin")
# ntc_cellchat |> saveRDS(
#   # TODO: change the path
#   here(
#     "processed-data/rds/spg_ccc",
#     "spd_neun_claudin_ntc_cellchat.rds"
#   )
# )

# rm(ntc_cellchat) # Clean Memory

# set.seed(20241107)
# scz_cellchat <- subset_N_cellchat(spe, .dx = "scz", var = "spd_neun_claudin")
# scz_cellchat |> saveRDS(
#   here(
#     "processed-data/rds/spg_ccc",
#     "spd_neun_claudin_scz_cellchat.rds"
#   )
# )

# ntc_cellchat <- readRDS(here(
#   "processed-data/rds/spg_ccc",
#   "spd_neun_claudin_ntc_cellchat.rds"
# ))
# scz_cellchat <- readRDS(
#   here(
#     "processed-data/rds/spg_ccc",
#     "spd_neun_claudin_scz_cellchat.rds"
#   )
# )

# object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)
# cellchat <- mergeCellChat(
#   object.list,
#   add.names = names(object.list)
# )

# saveRDS(
#   cellchat,
#   here(
#     "processed-data/rds/spg_ccc",
#     "spd_neun_claudin_merged_cellchat.rds"
#   )
# )







# Session Info ----
session_info()
