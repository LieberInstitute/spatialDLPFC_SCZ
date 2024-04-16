# Load packages ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(sessioninfo)
})


# Load Data ----
## Load Spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

## Load outlier keys ----
tot_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

# Keep non-outlier spots
spe <- spe[, spe$in_tissue == TRUE]
spe$all_outlier <- tot_outlier_df[spe$key, "all_outlier"]
spe <- spe[, spe$all_outlier == "neither"]

# Examine Genes ----
## (Deprecated) Genes with 0 count ----
# gene_sum <- rowSums(counts(spe))
# sum(gene_sum == 0)
# # [1] 3956
# mean(gene_sum == 0) # High proportion of non expressing gene
# # [1] 0.1080845
# gene_sum_0_index <- which(gene_sum != 0)
# spe <- spe[gene_sum != 0, ]
# stopifnot(all(rowSums(counts(spe)) != 0))


## Genes with 0 variance ----
# Note: Genes with 0 expression across samples is a subset
gene_var <- counts(spe) |> MatrixGenerics::rowVars()

sum(gene_var == 0)
# [1] 4061

gene_const_var <- which(gene_var == 0)
gene_exclude_df <- rowData(spe)[names(gene_const_var), ]

saveRDS(
  gene_exclude_df,
  here(
    "processed-data/rds/02_visium_qc",
    "excluded_genes_df.rds"
  )
)


# Session info ----
session_info()
