# Load library ----
suppressPackageStartupMessages({
  library(fgsea)
  # library(ggplot2)
  library(msigdbr) # Read in BrainGMT
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load data -----
## Load Pairwise test res ----
.spd <- "PRECAST_07"
layer_res <- readRDS(
  here(
    # TODO: need organize
    "processed-data", "rds", "layer_enrich_test",
    paste0("test_pairwise_", .spd, ".rds")
  )
)

## Load data base BrainGMT-----
BrainGMT <- gmtPathways(
  here(
    "code/analysis/SPG_analysis",
    "BrainGMTv2_HumanOrthologs.gmt.txt"
  )
)

# Annoate DEGs -----

## fgsea analysis ----
### BrainGMT ----
ordered_gene_vector_braingmt <- layer_res |>
  select(gene, ensembl, `t_stat_spd01-spd04`) |>
  arrange(desc(`t_stat_spd01-spd04`)) |>
  # remove genes whose have multiple ensemble
  filter(!duplicated(layer_res$gene, fromLast = TRUE)) |>
  select(-ensembl) |>
  deframe()

# error prevention
stopifnot(anyDuplicated(names(ordered_gene_vector_braingmt)) == 0)


set.seed(20250123)
braingmt_res <- fgsea(BrainGMT, ordered_gene_vector_braingmt,
  minSize = 10, maxSize = 1000,
  nPermSimple = 10000
)

braingmt_res$leadingEdge <- vapply(braingmt_res$leadingEdge,
  paste,
  collapse = ",", character(1L)
)

braingmt_res |>
  arrange(padj) |>
  filter(padj < 0.05) |>
  mutate(regulation = case_when(
    NES > 0 ~ "Enriched in SpD01 (L6/WM)",
    NES <= 0 ~ "Enriched in  SpD04 (WM)"
  )) |>
  write_csv(
    here(
      "processed-data/spd_marker_gene",
      "spd01-04-deg_brainGMT.csv"
    )
  )


# Session Info ----
session_info()
