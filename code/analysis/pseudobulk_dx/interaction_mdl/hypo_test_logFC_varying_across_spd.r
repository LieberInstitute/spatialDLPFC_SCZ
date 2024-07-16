# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
  library(ggrepel)
})

# Load interaction fit ----
fit <- readRDS(
  here(
    "processed-data/PB_dx_genes",
    "test_inter_PRECAST_07_20240627.rds"
  )
)

# Define inter terms ----
int_terms <- c(
  "dxscz:spdspd02", "dxscz:spdspd03",
  "dxscz:spdspd04", "dxscz:spdspd05", "dxscz:spdspd06",
  "dxscz:spdspd07"
)

# Hypo testing interaction effect ----
int_df <- topTable(
  fit,
  coef = int_terms, num = Inf
) |>
  rownames_to_column(var = "ensembl")

## Create gene names ---
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
) |> select(
  ensembl, gene
)

int_df <-  int_df |>
  left_join(
    gene_df,
    by = c("ensembl")
  )

## Save test results ----
saveRDS(
  int_df,
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_hypo_test_vary_logFC_across_spd.rds"
  )
)

# Session Info ----
sessioninfo::session_info()
