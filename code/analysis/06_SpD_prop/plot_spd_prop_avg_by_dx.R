# Average spatial domain (SpD) proportion: NTC vs SCZ
#
# For each sample, compute the proportion of spots assigned to each PRECAST_07
# SpD. Then average these per-sample proportions within each diagnosis group
# (each donor weighted equally), and visualize:
#   A) Stacked bar of the average SpD composition for NTC vs SCZ
#   B) Per-SpD comparison: box plot of per-sample proportions,
#      annotated with Wilcoxon rank-sum p-values (BH-adjusted across SpDs)

# Load library -----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(patchwork)
  library(sessioninfo)
})

# Load data ----
spe <- readRDS(
  here(
    "processed-data/rds",
    "03_visium_spatial_clustering",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

# Load PRECAST label ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds",
    "03_visium_spatial_clustering",
    "PRECAST_label_df_semi_sup_k_2-16.rds"
  )
)

## Merge PRECAST df ----
spe <- spe[, spe$key %in% PRECAST_df$key]
col_data_df <- PRECAST_df |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

# Load PRECAST07 annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  ),
  show_col_types = FALSE
) |>
  mutate(anno_lab = paste0(label, " (", spd, ")")) |>
  # laminar order: WM (bottom) -> L1/M (top)
  arrange(order)

# Calculate per-sample SpD proportion ----
spd_prop_df <- col_data_df |>
  transmute(
    sample_label,
    DX = toupper(dx),
    PRECAST_07
  ) |>
  dplyr::count(sample_label, DX, PRECAST_07, name = "n_spots") |>
  # Fill in SpDs absent from a sample with 0 spots, otherwise the group
  # averages are biased upward for those SpDs
  # (spd04/WM is absent in Br5182_NTC, Br5367_NTC, Br5436_NTC)
  complete(
    nesting(sample_label, DX), PRECAST_07,
    fill = list(n_spots = 0)
  ) |>
  group_by(sample_label) |>
  mutate(
    total_spots = sum(n_spots),
    proportion = n_spots / total_spots
  ) |>
  ungroup()

# error prevention
stopifnot(
  # every sample has all 7 SpDs after filling zeros
  all(dplyr::count(spd_prop_df, sample_label)$n == 7),
  # proportions sum to 1 with rounding error
  isTRUE(all.equal(
    spd_prop_df |>
      group_by(sample_label) |>
      summarize(s = sum(proportion)) |>
      pull(s),
    rep(1, n_distinct(spd_prop_df$sample_label)),
    tolerance = 1e-8
  ))
)

# Calculate average proportion per DX ----
spd_avg_df <- spd_prop_df |>
  group_by(DX, PRECAST_07) |>
  summarize(
    n_samples = n(),
    mean_prop = mean(proportion),
    sd_prop = sd(proportion),
    se_prop = sd_prop / sqrt(n_samples),
    ci_lo = mean_prop - qt(0.975, n_samples - 1) * se_prop,
    ci_hi = mean_prop + qt(0.975, n_samples - 1) * se_prop,
    median_prop = median(proportion),
    # pooled proportion (all spots of the group together), for reference
    pooled_prop = sum(n_spots) / sum(total_spots),
    .groups = "drop"
  )

# Average per-sample proportions sum to 1 within each DX
stopifnot(
  isTRUE(all.equal(
    spd_avg_df |>
      group_by(DX) |>
      summarize(s = sum(mean_prop)) |>
      pull(s),
    c(1, 1),
    tolerance = 1e-8
  ))
)

## Test difference between DX per SpD ----
spd_test_df <- spd_prop_df |>
  group_by(PRECAST_07) |>
  summarize(
    p_wilcox = wilcox.test(
      proportion[DX == "SCZ"], proportion[DX == "NTC"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) |>
  mutate(p_adj_BH = p.adjust(p_wilcox, method = "BH"))

## Combine summary table ----
spd_summary_df <- spd_avg_df |>
  pivot_wider(
    id_cols = PRECAST_07,
    names_from = DX,
    values_from = c(
      n_samples, mean_prop, sd_prop, se_prop, median_prop, pooled_prop
    )
  ) |>
  mutate(diff_mean_SCZ_minus_NTC = mean_prop_SCZ - mean_prop_NTC) |>
  left_join(spd_test_df, by = "PRECAST_07") |>
  left_join(
    spd_anno_df |> select(PRECAST_07 = spd, label, order),
    by = "PRECAST_07"
  ) |>
  arrange(order) |>
  relocate(label, .after = PRECAST_07) |>
  select(-order)

print(spd_summary_df, width = Inf)

write_csv(
  spd_summary_df,
  here(
    "code/analysis/06_SpD_prop",
    "spd_prop_avg_by_dx.csv"
  )
)

# Make plot ----
## Shared aesthetics ----
# Same SpD palette as plot_spd_proportion_ordered.R
spd_colors <- set_names(
  Polychrome::palette36.colors(7)[seq.int(7)],
  sprintf("spd%02d", 1:7)
)
dx_colors <- c("NTC" = "blue", "SCZ" = "red")

# Map spd -> annotated label, in laminar order (WM first)
spd_lab <- set_names(spd_anno_df$anno_lab, spd_anno_df$spd)

## A: Stacked bar of average composition ----
bar_df <- spd_avg_df |>
  mutate(
    # first factor level is stacked on top: L1/M on top, WM at bottom
    PRECAST_07 = factor(PRECAST_07, levels = rev(spd_anno_df$spd)),
    pct_lab = scales::percent(mean_prop, accuracy = 0.1),
    # white text on dark fills, black text on light fills
    txt_col = if_else(
      colSums(col2rgb(spd_colors[as.character(PRECAST_07)]) *
        c(0.299, 0.587, 0.114)) < 128,
      "white", "black"
    )
  )

n_per_dx <- spd_avg_df |> distinct(DX, n_samples)
dx_lab <- set_names(
  paste0(n_per_dx$DX, "\n(n = ", n_per_dx$n_samples, ")"),
  n_per_dx$DX
)

bar_p <- bar_df |>
  ggplot(aes(x = DX, y = mean_prop, fill = PRECAST_07)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = pct_lab, color = txt_col, group = PRECAST_07),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_color_identity() +
  scale_fill_manual(
    name = "Spatial Domain",
    values = spd_colors,
    breaks = rev(spd_anno_df$spd), # legend top-to-bottom matches the bars
    labels = spd_lab
  ) +
  scale_x_discrete(labels = dx_lab) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Average percentage of spots",
    title = "Average SpD composition"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "left"
  )

## B: Per-SpD NTC vs SCZ comparison ----
p_lab_df <- spd_test_df |>
  mutate(
    p_lab = paste0(
      "p = ", signif(p_wilcox, 2),
      "\nadj. p = ", signif(p_adj_BH, 2)
    )
  )

facet_levels <- spd_lab[rev(spd_anno_df$spd)] # L1/M first, WM last

box_p <- spd_prop_df |>
  mutate(spd_lab = factor(spd_lab[PRECAST_07], levels = facet_levels)) |>
  ggplot(aes(x = DX, y = proportion, color = DX)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA # outliers are already shown by the jittered points
  ) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1.2) +
  geom_text(
    data = p_lab_df |>
      mutate(spd_lab = factor(spd_lab[PRECAST_07], levels = facet_levels)),
    aes(x = 1.5, y = Inf, label = p_lab),
    inherit.aes = FALSE,
    vjust = 1.2, size = 2.5
  ) +
  facet_wrap(~spd_lab, nrow = 1) +
  scale_color_manual(values = dx_colors, guide = "none") +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0.02, 0.25)) # headroom for p-value label
  ) +
  labs(
    x = NULL,
    y = "Percentage of spots per sample",
    title = "Per-SpD percentage per sample"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

## Combine ----
bar_p + box_p +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(tag_levels = "A")

# Save plot -----
ggsave(
  here(
    "plots/03_visium_spatial_clustering",
    "bar_dot_plot_spd_prop_avg_by_dx.pdf"
  ),
  width = 14, height = 5
)

# Session Info ----
session_info()
