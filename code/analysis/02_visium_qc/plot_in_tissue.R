# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(scater)
library(spatialLIBD)
library(tidyverse)
library(sessioninfo)

# File Paths --------------------------------------------------------------
path_raw_spe <- here(
  "processed-data/rds",
  "01_build_spe",
  "test_raw_spe_w_spg_N63.rds")

# QC Code -----------------------------------------------------------------
raw_spe <- readRDS(
  path_raw_spe
)

spe <- raw_spe

# Test only
# spe <- raw_spe[, raw_spe$sample_id %in% unique(raw_spe$sample_id)[1:3]]

spe$in_tissue_fac <- factor(spe$in_tissue, levels = c(FALSE, TRUE), labels = c("Out", "In"))
levels(spe$in_tissue_fac)


vis_grid_clus(
  spe, 
  clustervar = "in_tissue_fac",
  pdf_file = here("plots/02_visium_qc/in_tissue_grid.pdf"),
  # alpha = 0.5,
  color = c("1" = "#d8b365", "2" = "#f5f5f5"),
  width = 7.5,
  height = 9,
  alpha = 0.9,
  point_size = 0.6 
)

# TODO:
# Problem:
# If we decide to include in the supplementary figure. we'll have to optimize the
# plotting to be compatible with small visuals.


# Session Info ------
session_info()