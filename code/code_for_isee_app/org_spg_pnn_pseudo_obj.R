# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(here)
  library(sessioninfo)
})

# Load Pseudobulked data ----
raw_sce <- readRDS(
  here(
    "processed-data/rds/PB_dx_spg",
    "pseudo_pnn_pos_donor_spd.rds"
  )
)

## Fitler out unused domains ----
ret_sce <- raw_sce[
  ,
  raw_sce$registration_variable %in%
    sprintf("spd%02d", c(2, 3, 5, 6, 7))
]

# Adjust data ----
## Change SpD07-L1 to SpD07-L1/M ----
ret_sce$spd_label <- factor(
  ret_sce$registration_variable,
  levels = c(
    "spd07",
    "spd06",
    "spd02",
    "spd05",
    "spd03"
  ),
  labels = c(
    "SpD07-L1/M",
    "SpD06-L2/3",
    "SpD02-L3/4",
    "SpD05-L5",
    "SpD03-L6"
  )
)

## Keep only relavent columns ----
col_df <- colData(ret_sce) |> as.data.frame()

subset_col_df <- col_df |>
  select(
    sample_id,
    DX,
    sex,
    age,
    rin,
    spd_label
  )

colData(ret_sce) <- DataFrame(subset_col_df)

ret_sce

## Remove outdated reduced dims ----
reducedDim(ret_sce, type = "10x_pca") <- NULL
reducedDim(ret_sce, type = "10x_tsne") <- NULL
reducedDim(ret_sce, type = "10x_umap") <- NULL
reducedDimNames(ret_sce)

# Save object ----
saveRDS(
  ret_sce,
  here(
    "processed-data/rds/99_pb_data_for_isee",
    "sce_pseudo_SPG_PNN_donor_spd.rds"
  )
)

# Session info ----
sessioninfo::session_info()