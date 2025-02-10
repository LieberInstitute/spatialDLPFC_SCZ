# Load library -----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
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
  summarize(n_spots = n()) |>
  ungroup() |>
  left_join(total_spots_per_sample, by = "sample_label") |>
  mutate(proportion = n_spots / total_spots) |>
  select(-total_spots)

# Explore samples that have missing spd
samples_w_missing_spd <- n_spots_per_spd |>
  group_by(sample_label) |>
  summarize(n_spd = n()) |>
  filter(n_spd != 7)

# sample_label n_spd
#   <chr>        <int>
# 1 Br5182_NTC       6
# 2 Br5367_NTC       6
# 3 Br5436_NTC       6


# make plot ----
## Calculate order ----
# calculate the sum of WM (spd01 and spd04)
tmp_long_df <- n_spots_per_spd |> pivot_wider(
  id_cols = "sample_label",
  names_from = PRECAST_07, values_from = proportion
)

# spd04 (WM) are missing for the following three samples
tmp_long_df |>
  filter(if_any(everything(), is.na))
# A tibble: 3 × 8
#   sample_label   spd01 spd02  spd03 spd04 spd05 spd06 spd07
#   <chr>          <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 Br5182_NTC   0.0385  0.288 0.224     NA 0.156 0.162 0.131
# 2 Br5367_NTC   0.0375  0.361 0.0744    NA 0.128 0.214 0.186
# 3 Br5436_NTC   0.00801 0.413 0.0691    NA 0.209 0.148 0.153

# error prevention
# proportion sums to 1 with rounding error
stopifnot(
  all(
    all.equal(
      rowSums(
        # remove sample_label column
        tmp_long_df[, -1] |>
          data.matrix(),
        na.rm = TRUE
      ),
      rep(1, nrow(tmp_long_df)),
      tolerance = 1e-8
    )
  )
)

# vec of sample_label ordered by WM proportion
sample_WM_ordered <- tmp_long_df |>
  replace_na(list(spd04 = 0)) |>
  mutate(
    WM_sum = rowSums(
      across(c(spd01, spd04)) # ,
      # na.rm = TRUE
    )
  ) |>
  arrange(WM_sum) |>
  pull(sample_label)

## Make plots ----
ret_p <- n_spots_per_spd |>
  mutate(PRECAST_07 = factor(PRECAST_07, levels = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4)))) |>
  ggplot() +
  geom_bar(
    aes(x = sample_label, y = proportion, fill = PRECAST_07),
    stat = "identity",
    position = "fill"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(
    x = "Sample Label",
    y = "Proportion",
    fill = "PRECAST_07"
  ) +
  scale_x_discrete(
    # change the order of the x-axis
    limits = sample_WM_ordered
  )



# Change color palette and spd ordering
# Load PRECAST07 annotation
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

spd_order <- order(spd_anno_df$anno_lab)

bar_p <- ret_p +
  scale_fill_manual(
    name = "Spatial Domain",
    # change fill palette
    values = set_names(
      Polychrome::palette36.colors(7)[seq.int(7)],
      unique(spe$PRECAST_07) |> sort()
    ),
    # change order of fill
    breaks = spd_anno_df$spd[spd_order],
    labels = spd_anno_df$anno_lab[spd_order]
  ) +
  theme_minimal()

## Create symbol for sample id ----
demo_df <- metadata(spe)$dx_df |>
  transmute(
    DX = toupper(dx),
    sample_label = paste0(
      subject, "_", toupper(dx)
    ),
    sex,
    h_pos = 1,
  )


dx_p <- demo_df |>
  ggplot() +
  geom_point(
    aes(
      x = sample_label,
      y = h_pos,
      color = DX,
      shape = sex
    ),
    size = 3
  ) +
  scale_x_discrete(
    # change the order of the x-axis
    limits = sample_WM_ordered
  ) +
  scale_color_manual(
    values = c("NTC" = "blue", "SCZ" = "red"),
    guide = "none"
  )
  # ) #+
# theme_void()

tmp_p <- dx_p +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


library(cowplot)

plot_grid(
  # Adjust plot
  bar_p +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    ),
  tmp_p +
    # theme_void() +
    theme(
      # legend.position = "none",
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA)
    ),
  align = "v",
  ncol = 1,
  rel_heights = c(0.8, 0.2)
)


# Save plot -----
ggsave(
  here(
    "plots/03_visium_spatial_clustering",
    "bar_plot_spd_prop_per_sample.pdf"
  ),
  width = 11, height = 5,
)

# Session Info ----
sessioninfo()
