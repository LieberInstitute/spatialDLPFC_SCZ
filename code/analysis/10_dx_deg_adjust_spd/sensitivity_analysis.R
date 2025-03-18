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
    "processed-data/rds/10_dx_deg_adjust_spd/preliminary",
    "test_PRECAST_07.csv"
  )
)

## Load new DEG res ----
RIN_adj_res <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "tmp_rin_adj_PRECAST07.csv"
  )
)

## Merged data sets ----
nrow(prelim_res) == nrow(RIN_adj_res)

merged_res <- inner_join(
  prelim_res,
  RIN_adj_res,
  by = c("ensembl"),
  suffix = c("_prelim", "_rin_adj")
)


prelim_res |> filter(gene == "GBP2")

# histogram of p-value
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
# [1] 0.6997641

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
# [1] 0.7198662

# Session Info ----
sessioninfo::session_info()
