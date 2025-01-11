# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(scater)
  library(sessioninfo)
})

# Load Data ----
## Load Spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "raw_spe_wo_SPG_N63.rds"
  )
)

spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$brnum <- metadata(spe)$dx_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sample_id <- paste0(
  spe$brnum
)


## Load outlier keys ----

tot_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

## Merge data together ----
spe$all_outlier <- tot_outlier_df[spe$key, "all_outlier"]
spe$discard <- tot_outlier_df[spe$key, "remove"]


## Remove out_tissue_spots ----
spe <- spe[, spe$in_tissue == TRUE]

# Create data frame ----
##
sample_lvl_summary_df <- colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(
    n_total = n(),
    n_discard = sum(discard),
    prop_discard = n_discard / n_total
  ) |>
  right_join(
    metadata(spe)$dx_df,
    by = join_by(sample_id == subject)
  )

## Outlier samples ---
# Expectation: samples with regional artifacts ----

sample_lvl_summary_df |>
  filter(prop_discard > 0.05) |>
  select(sample_id, sample_id.y)
# A tibble: 2 × 2
# sample_id sample_id.y
# <chr>     <chr>
# 1 Br1113    V13M06-279_A1
# 2 Br6437    V13M06-282_B1

# Box plot of Dx groups----
pdf(
  here("plots/02_visium_qc/", "num_prop_discard_spots.pdf")
)
## Num of spots ----
sample_lvl_summary_df |>
  ggplot(aes(x = dx, y = n_discard)) +
  geom_boxplot(aes(color = dx)) +
  geom_jitter(size = 0.5, alpha = 0.3) +
  theme_light() +
  ylab("# of Disarded Spots") +
  xlab("") +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  )

## Prop of spots ----
# NOTE: out of in-tissue spots
sample_lvl_summary_df |>
  ggplot(aes(x = dx, y = prop_discard)) +
  geom_boxplot(aes(color = dx)) +
  geom_jitter(size = 0.5, alpha = 0.3) +
  theme_light() +
  ylab("Prop of Disarded Spots") +
  xlab("") +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  )

dev.off()

# Cumulative curves for QC metrics ----
discard_spe <- spe[, spe$discard == TRUE]

pdf(
  here("plots/02_visium_qc/", "distribution_discard_spots.pdf")
)

## UMI ----
colData(discard_spe) |>
  data.frame() |>
  ggplot() +
  stat_ecdf(
    aes(x = sum_umi, group = sample_id, color = dx),
    alpha = 0.5,
    geom = "step"
  ) +
  coord_cartesian(xlim = c(1, 2000)) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Cumulative Distribution of UMI Counts",
    x = "Log10(UMI Counts)",
    y = "Cumulative Probability",
    color = "Diagnosis"
  ) +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

## Unique genes ----
colData(discard_spe) |>
  data.frame() |>
  ggplot() +
  stat_ecdf(
    aes(x = sum_gene, group = sample_id, color = dx),
    alpha = 0.5,
    geom = "step"
  ) +
  # scale_x_log10() +
  coord_cartesian(xlim = c(1, 4000)) +
  theme_minimal() +
  labs(
    title = "Cumulative Distribution of Unique Genes",
    x = "Log10(Unique Genes)",
    y = "Cumulative Probability",
    color = "Diagnosis"
  ) +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

## Mito Ratio ----
colData(discard_spe) |>
  data.frame() |>
  ggplot() +
  stat_ecdf(
    aes(x = expr_chrM_ratio, group = sample_id, color = dx),
    alpha = 0.5,
    geom = "step"
  ) +
  theme_minimal() +
  labs(
    title = "Cumulative Distribution of Mitochondrial Ratio",
    x = "Mitochondrial Ratio",
    y = "Cumulative Probability",
    color = "Diagnosis"
  ) +
  scale_color_manual(
    values = c(
      "ntc" = "blue",
      "scz" = "red"
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )
dev.off()

# Session Info ----
sessioninfo::session_info()
