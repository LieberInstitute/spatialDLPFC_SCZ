# Load Library -----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(gtsummary) # Calculate quick summary stat for spd prop
  library(sessioninfo)
})

# Load data ----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

# Load PRECAST label ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Merge PRECAST df ----
precast_vars <- grep(
  "^PRECAST_", colnames(PRECAST_df),
  value = TRUE
)

spe <- spe[, spe$key %in% PRECAST_df$key]
col_data_df <- PRECAST_df |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )
rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# Calculate composition df ----
total_spots_per_sample <- col_data_df |>
  group_by(sample_label) |>
  summarize(total_spots = n())

n_spots_per_spd <- col_data_df |>
  group_by(
    sample_label,
    PRECAST_07
  ) |>
  summarize(
    n_spots = n(), dx = unique(dx)
  ) |>
  ungroup() |>
  left_join(total_spots_per_sample, by = "sample_label") |>
  mutate(
    proportion = n_spots / total_spots,
    DX = toupper(dx)
  ) |>
  select(-total_spots)

# Explore samples that have missing spd
# samples_w_missing_spd <- n_spots_per_spd |>
#   group_by(sample_label) |>
#   summarize(n_spd = n()) |>
#   filter(n_spd != 7)

# sample_label n_spd
#   <chr>        <int>
# 1 Br5182_NTC       6
# 2 Br5367_NTC       6
# 3 Br5436_NTC       6


## Calculate order ----
# # calculate the sum of WM (spd01 and spd04)
tmp_long_df <- n_spots_per_spd |> pivot_wider(
  id_cols = c("sample_label", "DX"),
  names_from = PRECAST_07, values_from = proportion
)

# # spd04 (WM) are missing for the following three samples
# tmp_long_df |>
#   filter(if_any(everything(), is.na))
# A tibble: 3 × 8
#   sample_label   spd01 spd02  spd03 spd04 spd05 spd06 spd07
#   <chr>          <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Br5182_NTC   0.0385  0.288 0.224     NA 0.156 0.162 0.131
# 2 Br5367_NTC   0.0375  0.361 0.0744    NA 0.128 0.214 0.186
# 3 Br5436_NTC   0.00801 0.413 0.0691    NA 0.209 0.148 0.153

tmp_long_df <- tmp_long_df |>
  replace_na(list(spd04 = 0)) |>
  mutate(
    WM_sum = rowSums(
      across(c(spd01, spd04)) # ,
      # na.rm = TRUE
    )
  )

# Make plot ----
## Make spd-facet plot ----
# Find data in long form
n_spots_per_spd |>
  ggplot() +
  geom_boxplot(aes(x = DX, y = proportion, fill = DX)) +
  facet_wrap(~PRECAST_07) +
  theme_minimal() +
  scale_fill_manual(
    values = c("NTC" = "blue", "SCZ" = "red")
  )
# TODO: Order the facet based on spatial domain


# Wilcoxon tests (not included in the manuscript)
tmp_long_df |>
  tbl_summary(include = c(starts_with("spd"), WM_sum), by = DX) |>
  add_p(test = list(all_continuous() ~ "wilcox.test"))

# Characteristic	NTC
# N = 31	SCZ
# N = 32	p-value
# spd01	0.08 (0.05, 0.09)	0.06 (0.05, 0.08)	0.4
# spd02	0.35 (0.25, 0.38)	0.33 (0.26, 0.41)	>0.9
# spd03	0.14 (0.09, 0.18)	0.14 (0.08, 0.17)	0.9
# spd04	0.04 (0.02, 0.10)	0.06 (0.02, 0.11)	0.6
# spd05	0.11 (0.09, 0.13)	0.12 (0.09, 0.14)	0.6
# spd06	0.14 (0.09, 0.16)	0.13 (0.10, 0.15)	0.4
# spd07	0.14 (0.09, 0.17)	0.12 (0.08, 0.17)	0.6
# WM_sum	0.12 (0.09, 0.19)	0.13 (0.08, 0.19)	0.8
# 1 Median (Q1, Q3)
# 2 Wilcoxon rank sum exact test; Wilcoxon rank sum test


## Make chart for marginal spd prop
n_spots_per_spd |>
  ggplot() +
  geom_boxplot(aes(x = PRECAST_07, y = proportion)) +
  theme_minimal()
# TODO: add color to spatial domain
# TODO: Adjust the order of the spatial domains.

# calculate the statistics of prop of SpD

tmp_long_df |>
  tbl_summary(include = c(starts_with("spd"), WM_sum))

# Characteristic	N = 63
# spd01	0.07 (0.05, 0.09)
# spd02	0.34 (0.25, 0.41)
# spd03	0.14 (0.08, 0.17)
# spd04	0.05 (0.02, 0.11)
# spd05	0.12 (0.09, 0.14)
# spd06	0.13 (0.09, 0.16)
# spd07	0.13 (0.08, 0.17)
# WM_sum(spd01+spd04)	0.12 (0.09, 0.19)
# 1 Median (Q1, Q3)


# Session Info ----
session_info()
