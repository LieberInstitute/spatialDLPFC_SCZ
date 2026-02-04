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

## Significant at FDR < 0.10 ----
test_res |>
  filter(adj.P.Val < 0.10) |> nrow()

test_res |>
  filter(adj.P.Val < 0.10) |>
  select(gene_name, logFC, P.Value, adj.P.Val) |>
  View()



test_res |>
  filter(gene_name %in% c("FKBP5", "XRRA1")) |> View()



# Session Info ----
session_info()
