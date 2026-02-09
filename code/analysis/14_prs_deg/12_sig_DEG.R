# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(sessioninfo)
  library(here)
})

# Load data ----
test_res <- read_csv(
  file = here(
    "processed-data/rds/14_prs_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)

# Some statistics ----
## Significant at FDR < 0.10 ----
test_res |>
  filter(adj.P.Val < 0.10) |> nrow()

## Significant at nom p-value  ----
test_res |>
  filter(P.Value < 0.05) |> nrow()

test_res |>
#   filter(adj.P.Val < 0.10) |>
#   select(gene_name, logFC, P.Value, adj.P.Val) |>
#   View()

test_res |>
  filter(
    gene_name %in% c(
      "EIF2A", "COX7A2", "MCL1", "AIF1", "TREM2",
      "VEGFA", "A2M", "FOSB", "JUN", "SYT1", "KANSL1-AS1",
      "MAP2", "GAP43", "MAPK3", "MSI2", "R3HDM2", "THOC7",
       "C1QA", "C1QB", "C3", "CX3CR1", "TYROBP", "CD74")) |> 
    View()

# Load Dx-DEG results ----
dx_deg_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)
dx_deg_df |> dim()

dx_deg_df |>
  filter(fdr_scz < 0.10) |> nrow()


# Overlap between Dx-DEGs and PRS-DEGs ----
overlap_ensmbl <- intersect(
  dx_deg_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
  test_res |> filter(adj.P.Val < 0.10) |> pull(gene_id)
)

overlap_ensmbl |> length()

overlap_gene_name <- test_res |>
  filter(gene_id %in% overlap_ensmbl) |>
  pull(gene_name)

stopifnot(all(c(
  "AIF1", "TREM2", "C1QA", "C1QB",
   "C3", "CX3CR1", "TYROBP", "CD74") %in%
    overlap_gene_name))


# Session Info ----
session_info()
