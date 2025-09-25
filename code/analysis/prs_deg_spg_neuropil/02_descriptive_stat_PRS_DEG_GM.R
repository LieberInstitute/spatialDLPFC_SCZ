# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
})

# Load Neuropil PRS-DEG results ----
prs_deg_res <- read_csv(
  here(
      "processed-data/rds/14_prs_deg",
      "Neuropil_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
    )
)

# Descriptive statistics ----
## Histogram of p-values ----
hist(
  prs_deg_res$P.Value,
  breaks = 100,
  main = "Histogram of p-values",
  xlab = "p-value",
  ylab = "Frequency"
)

## Number of DEGs at diffferent p-cutoffs ----
### Nominal p < 0.05 ----
sum(prs_deg_res$P.Value < 0.05)
# [1] 354

### FDR < 0.1 ----
sum(prs_deg_res$adj.P.Val < 0.1)
# [1] 5

prs_deg_res |> filter(adj.P.Val < 0.1) |> pull(gene_name)
# [1] "MTRNR2L8"   "AC004556.3" "SURF1"      "MRPL23"    
# [5] "C2orf74" 



# Session info ----
session_info()