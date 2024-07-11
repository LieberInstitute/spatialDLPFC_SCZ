library(SpatialExperiment)
library(here)
library(tidyverse)

spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)


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

# Box plot ----
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
  geom_jitter(point = 0.5, alpha = 0.3) +
  theme_light() +
  ylab("# of Spots") +
  xlab("") +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  )
