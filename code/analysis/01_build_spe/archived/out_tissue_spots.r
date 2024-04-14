# Load Pacakges ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(scater)
  library(ggpubr)
  library(sessioninfo)
})

# Load Data ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63.rds"
  )
)

# mito% and sum gene relationship ----
pdf(here("plots/02_visium_qc/compare_in_out_tissue_mito_gene.pdf"))
## Overall ----
p_all_mito <- plotColData(
  # Plot in_tissue=FALSE spots last
  spe[, order(colData(spe)$in_tissue, decreasing = TRUE)],
  x = "sum_gene", y = "expr_chrM_ratio",
  colour_by = "in_tissue",
  point_size = 0.05
) + labs(title = "All samples")

p_all_umi <- plotColData(
  # Plot in_tissue=FALSE spots last
  spe[, order(colData(spe)$in_tissue, decreasing = TRUE)],
  x = "sum_gene", y = "sum_umi",
  colour_by = "in_tissue",
  point_size = 0.05
) + labs(title = "All samples")

ggpubr::ggarrange(p_all_mito, p_all_umi) |>
  print()



## Per-sample ----
for(.smp in unique(spe$sample_id)){
  sub_spe <- spe[, spe$sample_id == .smp]
  lay_order <- order(colData(sub_spe)$in_tissue, decreasing = TRUE)

  p_mito <- plotColData(
    sub_spe[, lay_order],
    x = "sum_gene", y = "expr_chrM_ratio",
    colour_by = "in_tissue"
  ) + labs(title = .smp)

  p_umi <- plotColData(
    sub_spe[, lay_order],
    x = "sum_gene", y = "sum_umi",
    colour_by = "in_tissue"
  ) + labs(title = .smp)

  p <- ggpubr::ggarrange(p_mito, p_umi)
  print(p)

}
dev.off()

# Session ----
session_info()
