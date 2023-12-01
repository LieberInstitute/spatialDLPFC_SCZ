
# Load Packages -----------------------------------------------------------
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)

# Load SPE ----------------------------------------------------------------
path_spe_after_spot_qc <- here::here(
  "processed-data", "rds",
  "01_build_spe/",
  "raw_spe_wo_SPG_N63.rds")

# Plot all samples --------------------------------------------------------

vis_grid_gene(
  spe,
  geneid = "sum_gene",
  spatial = FALSE,
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_gene.pdf"),
  sample_order = unique(spe$sample_id) |> sort(),
  point_size = 0.8)

vis_grid_gene(
  spe,
  geneid = "sum_umi",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_umi.pdf"),
  sample_order = unique(spe$sample_id) |> sort(),
  point_size = 0.8)

vis_grid_gene(
  spe,
  geneid = "expr_chrM",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_expr_chrM.pdf"),
  sample_order = unique(spe$sample_id) |> sort(),
  point_size = 0.8 )

vis_grid_gene(
  spe,
  geneid = "expr_chrM_ratio",
  pdf_file = here("plots/02_visium_qc/test_oddity",
                  "random_pattern_expr_chrM_ratio.pdf"),
  sample_order = unique(spe$sample_id) |> sort(),
  point_size = 0.8 )

# Problematic Samples Only ------------------------------------------------

prob_spe <- spe[, str_starts(spe$sample_id, "V12F14")]

vis_grid_gene(
  prob_spe,
  geneid = "sum_gene",
  spatial = FALSE,
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_gene-prob_only.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)

vis_grid_gene(
  prob_spe,
  geneid = "sum_umi",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_umi-prob_only.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)

vis_grid_gene(
  prob_spe,
  geneid = "expr_chrM",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_expr_chrM-prob_only.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)

vis_grid_gene(
  prob_spe,
  geneid = "expr_chrM_ratio",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_expr_chrM_ratio-prob_only.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)

# Proof showing good after log normalization ------------------------------




# Session Info ----
sessioninfo::session_info()
