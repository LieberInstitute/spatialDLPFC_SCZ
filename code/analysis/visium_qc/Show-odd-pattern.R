# Load Packages -----------------------------------------------------------
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(scuttle)
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

## Log lib size norm ----------------------------------------------------
prob_spe <- logNormCounts(prob_spe)
prob_spe <- logNormCounts(prob_spe, log = FALSE)


## Normalized Plot ---------------------------------------------------------
prob_spe$norm_sum_umi <- assay(prob_spe, "normcounts") |> colSums()
vis_grid_gene(
  prob_spe,
  geneid = "norm_sum_umi",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_umi-prob_only-norm.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)
# Note: this is not informative at all, because every spots will be normalized
#       to have the same UMI

## Lognormalized plot  -----------------------------------------------------
prob_spe$lognorm_sum_umi <- logcounts(prob_spe) |> colSums()
vis_grid_gene(
  prob_spe,
  geneid = "lognorm_sum_umi",
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_umi-prob_only-lognorm.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)

## Lognormalized plot (MOBP)  -----------------------------------------------------
vis_grid_gene(
  prob_spe,
  geneid = rowData(spe)$gene_id[which(rowData(spe)$gene_name == "MOBP")],
  pdf_file = here(
    "plots/02_visium_qc/test_oddity",
    "random_pattern_sum_umi-prob_only-lognorm-MOBP.pdf"),
  sample_order = unique(prob_spe$sample_id) |> sort(),
  point_size = 3)
# TODO: maybe need to change to something else.




# Session Info ----
sessioninfo::session_info()
