# Load Libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
})

# Load Data ----
## Load PRS DEGs results ----
test_res <- read_csv(
  file = here(
    "processed-data/rds/14_prs_deg",
    "PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)

# Plot histogram of p-values ----

hist(
  test_res$P.Value,
  breaks = 20, # bin width = 0.05
  main = "Histogram of -log10(p-value) for PRS DEGs",
  xlab = "-log10(p-value)",
  col = "lightblue",
  border = "black"
)


