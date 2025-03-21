# Load Library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load Data ----
## Load preliminary DEG res ----
prelim_res <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)

## Load new DEG res ----
RIN_adj_res <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_rin_adj_PRECAST07.csv"
  )
)
# Error prevention
stopifnot(nrow(prelim_res) == nrow(RIN_adj_res))

## Merged data sets ----
merged_res <- inner_join(
  prelim_res,
  RIN_adj_res,
  by = c("ensembl"),
  suffix = c("_prelim", "_rin_adj")
)

# histogram of p-value ----
## prelim_res ----
hist(prelim_res$p_value_scz, breaks = 100)
qqplot(prelim_res$p_value_scz, runif(n = nrow(prelim_res)))
abline(a = 0, b = 1, col = "red")

## RIN_adj_res ----
hist(RIN_adj_res$p_value_scz, breaks = 100)
qqplot(RIN_adj_res$p_value_scz, runif(n = nrow(RIN_adj_res)))
abline(a = 0, b = 1, col = "red")

# Sensitivity Analysis ----
## Scatter plot of LogFC ----
plot(
  merged_res$logFC_scz_prelim,
  merged_res$logFC_scz_rin_adj,
  xlab = "LogFC Preliminary",
  ylab = "LogFC RIN Adjusted",
  main = "LogFC Comparison"
)
abline(a = 0, b = 1, col = "red")

cor(
  merged_res$logFC_scz_prelim,
  merged_res$logFC_scz_rin_adj
)
# [1] 0.9578429

## Scatter plot of t-statistics ----
plot(
  merged_res$t_stat_scz_prelim,
  merged_res$t_stat_scz_rin_adj,
  xlab = "t-stat Preliminary",
  ylab = "t-stat RIN Adjusted",
  main = "t-stat Comparison"
)
abline(a = 0, b = 1, col = "red")

cor(
  merged_res$t_stat_scz_prelim,
  merged_res$t_stat_scz_rin_adj
)
# [1] 0.9563269


sum(merged_res$fdr_scz_rin_adj < 0.10)


## What are overlapping sig genes
merged_res |>
  mutate(
    prelim_sig = fdr_scz_prelim < 0.10,
    rin_sig = fdr_scz_rin_adj < 0.10) |>
  with(
    table(prelim_sig, rin_sig)
  ) # 81 genes

#             rin_sig
# prelim_sig FALSE  TRUE
#      FALSE 13763    42
#      TRUE     69   103

## RIN specific genes
merged_res |>
  filter(
    fdr_scz_prelim > 0.10,
    fdr_scz_rin_adj < 0.10
  ) |>
  pull(
    gene_rin_adj
  )

## Preliminary specific genes
merged_res |>
  filter(
    fdr_scz_prelim < 0.10,
    fdr_scz_rin_adj > 0.10
  ) |>
  pull(
    gene_prelim
  )


# Session Info ----
sessioninfo::session_info()
