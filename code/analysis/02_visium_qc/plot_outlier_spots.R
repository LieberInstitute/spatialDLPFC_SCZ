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

## Load local outlier keys ----
local_outlier_df <- readRDS(
  here(
    "processed-data/02_visium_qc",
    "local_outlier_df.rds"
  )
)




## Merge data together ----

tot_outlier_d <- full_join(
  outlier_df,
  local_outlier_df,
  by = "key"
)

tot_outlier_d <- tot_outlier_d |>
  rowwise() |>
  mutate(
    outlier =
      any(
        c(umi_lt_100, gene_lt_200, artifact)
      )
  ) |>
  ungroup() 
  
tot_outlier_d <- tot_outlier_d |>
  mutate(
    all_outlier = case_when(
      outlier & local_outliers ~ "Both",
      outlier & !local_outliers ~ "outlier",
      !outlier & local_outliers ~ "local",
      !outlier & !local_outliers ~ "neither"
    ) |> factor()
  ) |> column_to_rownames("key")

# Quick test
tot_outlier_d |> filter(local_outliers == TRUE) |> 
select(outlier, local_outliers, all_outlier)

spe$all_outlier <- tot_outlier_d[spe$key, "all_outlier"]

# Make spot plot ----
## Remove out-tissue spots ----
spe <- spe[, spe$in_tissue == TRUE]

vis_grid_clus(
  spe,
  clustervar = "all_outlier",
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
