# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(sessioninfo)
})

# Load Data ----
## Load Spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63.rds"
  )
)

## Load outlier keys ----
outlier_df <- readRDS(
  here(
    "processed-data/02_visium_qc",
    "outlier_df.rds"
  )
)

## Merge data together ----
outlier_keys <- outlier_df |>
  rowwise() |>
  mutate(
    outlier =
      any(
        c(umi_lt_100, gene_lt_200, artifact)
      )
  ) |>
  ungroup() |>
  filter(outlier == TRUE) |>
  pull(key)

spe$outlier <- factor(spe$key %in% outlier_keys)

## Remove out-tissue spots
spe <- spe[, spe$in_tissue == TRUE]

# Spot plot of outlier ----
vis_grid_clus(
  spe,
  clustervar = "outlier",
  spatial = FALSE,
  pdf_file = here(
    "plots/02_visium_qc/spot_plot_outliers.pdf"
  ),
  sample_order = unique(spe$sample_id) |> sort(),
  height = 1056,
  width = 816,
  point_size = 0.8
)

# Session Info ----
session_info()
