# Load library ----
library(here)
library(tidyverse)

# Load data ----
norm_prs_res <- read_csv(
  here(
    "processed-data/rds/14_prs_deg/norm_PRS_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)

prs_res <- read_csv(
  here(
    "processed-data/rds/14_prs_deg",
    "PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)

# Compare number of significant genes ----
norm_prs_sig_genes <- norm_prs_res |>
  filter(adj.P.Val < 0.10) |>
  nrow()

prs_sig_genes <- prs_res |>
  filter(adj.P.Val < 0.10) |>
  nrow()

norm_prs_sig_genes == prs_sig_genes
# TRUE


View(norm_prs_res)
View(prs_res)
