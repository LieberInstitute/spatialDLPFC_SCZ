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
    "pseudo_vasc_pos_donor_spd.rds"
  )
)

# Adjust data ----
ret_sce <- raw_sce

## Change SpD07-L1 to SpD07-L1/M ----
ret_sce$spd_label <- factor(
  ret_sce$PRECAST_07,
  levels = c(
    "spd07",
    "spd06",
    "spd02",
    "spd05",
    "spd03",
    "spd01",
    "spd04"
  ),
  labels = c(
    "SpD07-L1/M",
    "SpD06-L2/3",
    "SpD02-L3/4",
    "SpD05-L5",
    "SpD03-L6",
    "SpD01-WMtz",
    "SpD04-WM"
  )
)


## Keep only relavent columns ----
col_df <- colData(ret_sce) |> as.data.frame()

subset_col_df <- col_df |>
  mutate(
    sample_id = paste0(brnum, "_", toupper(dx)),
    DX = toupper(dx)
  ) |>
  left_join(
    metadata(raw_sce)$dx_df |> select(subject, rin = RIN),
    by = c("brnum" = "subject")
  ) |>
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

# Save the processed SCE object ----
saveRDS(
  ret_sce,
  here(
    "processed-data/rds/99_pb_data_for_isee",
    "sce_pseudo_SPG_vasc_donor_spd.rds"
  )
)

# Session info ----
session_info()
