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
    "raw_spe_wo_SPG_N63.rds"
  )
)

## Load outlier keys ----
outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "outlier_df.rds"
  )
)

## Load local outlier keys ----
local_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "local_outlier_df.rds"
  )
)




## Merge data together ----

tot_outlier_df <- full_join(
  outlier_df,
  local_outlier_df,
  by = "key"
)

tot_outlier_df <- tot_outlier_df |>
  rowwise() |>
  mutate(
    outlier =
      any(
        c(umi_lt_100, gene_lt_200, artifact)
      )
  ) |>
  ungroup()

tot_outlier_df <- tot_outlier_df |>
  mutate(
    all_outlier = case_when(
      outlier & local_outliers ~ "Both",
      outlier & !local_outliers ~ "outlier",
      !outlier & local_outliers ~ "local",
      !outlier & !local_outliers ~ "neither"
    ) |> factor(),
    remove = (out_tissue_spots) | (all_outlier != "neither")
  ) |>
  column_to_rownames("key")

saveRDS(
  tot_outlier_df,
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

# Quick test
tot_outlier_df |>
  filter(local_outliers == TRUE) |>
  select(outlier, local_outliers, all_outlier) |> 
  head()

spe$all_outlier <- tot_outlier_df[spe$key, "all_outlier"]

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
