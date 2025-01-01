# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
})

# Load data (outliers removed) ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)

# Characterizing statistics ----
## Total spots ----
ncol(spe)
# [1] 279806

## Total spots per dx ----
colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(n = n()) |>
  right_join(
    metadata(spe)$dx_df,
    by = "sample_id"
  ) |>
  group_by(dx) |>
  summarise(
    n = sum(n)
  )

#   dx         n
#   <chr>  <int>
# 1 ntc   136456
# 2 scz   143350

# Medians ----
colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(n = n()) |>
  right_join(
    metadata(spe)$dx_df,
    by = "sample_id"
  ) |>
  group_by(dx) |>
  summarize(
    median = median(n)
  )

#   dx    median
#   <chr>  <dbl>
# 1 ntc     4556
# 2 scz     4590


# Plot ----
## Box plot ----
colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(n = n()) |>
  right_join(
    metadata(spe)$dx_df,
    by = "sample_id"
  ) |>
  ggplot(aes(x = dx, y = n)) +
  geom_boxplot(aes(color = dx)) +
  geom_jitter(size = 0.5, alpha = 0.3) +
  theme_light() +
  ylab("# of Spots") +
  xlab("") +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  )

# Session Info ----
session_info()
