# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(here)
  library(ggpubr)
  library(sessioninfo)
})

# Load Pseudobulked data ----
raw_sce <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_PRECAST07_donor_spd.rds"
  )
)

# Adjust data ----
ret_sec <- raw_sce
## Change SpD07-L1 to SpD07-L1/M ----
levels(ret_sec$spd_label)[levels(ret_sec$spd_label) == "SpD07-L1"] <- "SpD07-L1/M"
levels(ret_sec$spd_label)

## Keep only relavent columns ----
col_df <- colData(raw_sce) |> as.data.frame()

subset_col_df <- col_df |>
  select(
    sample_id = sample_label,
    DX,
    sex,
    age,
    rin,
    spd_label
  )

colData(ret_sec) <- DataFrame(subset_col_df)

ret_sec

## Remove outdated reduced dims ----
reducedDim(ret_sec, type = "10x_pca") <- NULL
reducedDim(ret_sec, type = "10x_tsne") <- NULL
reducedDim(ret_sec, type = "10x_umap") <- NULL
reducedDimNames(ret_sec)

# Check reduced dim for the pseudobulked data ----
set.seed(20250725)
# Calculate PCA
ret_sec <- runPCA(ret_sec, exprs_values = "logcounts")
# Plot PCA using scater
plotPCA(ret_sec, colour_by = "spd_label") +
  ggtitle("PCA of Pseudobulked Data")

# Calculate UMAP
set.seed(20250725)
ret_sec <- runUMAP(ret_sec, dimred = "PCA")

ggarrange(
  plotReducedDim(
    ret_sec,
    dimred = "UMAP",
    colour_by = "spd_label"
  ) +
    ggtitle("UMAP of Pseudobulked Data"),
  plotReducedDim(
    ret_sec,
    dimred = "UMAP",
    colour_by = "sex"
  ) +
    ggtitle("UMAP of Pseudobulked Data"),
  ncol = 2
)

# Calculate tSNE
set.seed(20250725)
ret_sec <- runTSNE(ret_sec)

ggarrange(
  plotReducedDim(
    ret_sec,
    dimred = "TSNE",
    colour_by = "spd_label"
  ) +
    ggtitle("TSNE of Pseudobulked Data"),
  plotReducedDim(
    ret_sec,
    dimred = "TSNE",
    colour_by = "sex"
  ) +
    ggtitle("TSNE of Pseudobulked Data"),
  ncol = 2
)

# Save object ----
saveRDS(
  ret_sec,
  here(
    "processed-data/rds/99_pb_data_for_isee",
    "sce_pseudo_PRECAST07_donor_spd_for_isee.rds"
  )
)

# Session info ----
session_info()
