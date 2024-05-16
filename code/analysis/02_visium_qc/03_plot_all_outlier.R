# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
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

## Fetch demo info
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
  spe$brnum, "_", spe$dx
)

## Load outlier keys ----
outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "outlier_df.rds"
  )
)

## Load local outlier keys ----
local_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "local_outlier_df.rds"
  )
)




## Merge data together ----
tot_outlier_df <- full_join(
  outlier_df,
  local_outlier_df,
  by = "key"
)

tot_outlier_df <- tot_outlier_df |>
  rowwise() |>
  mutate(
    outlier =
      any(
        c(umi_lt_100, gene_lt_200, artifact)
      )
  ) |>
  ungroup()

tot_outlier_df <- tot_outlier_df |>
  mutate(
    all_outlier = case_when(
      outlier & local_outliers ~ "Both",
      outlier & !local_outliers ~ "outlier",
      !outlier & local_outliers ~ "local",
      !outlier & !local_outliers ~ "neither"
    ) |> factor(),
    remove = (out_tissue_spots) | (all_outlier != "neither")
  ) |>
  column_to_rownames("key")


saveRDS(
  tot_outlier_df,
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

# Quick test
# tot_outlier_df |>
#   filter(local_outliers == TRUE) |>
#   select(outlier, local_outliers, all_outlier) |>
#   head()

spe$all_outlier <- tot_outlier_df[spe$key, "all_outlier"]
spe$discard <- tot_outlier_df[spe$key, "remove"]

# Make spot plot ----
spe <- spe[, spe$in_tissue == TRUE]

# Fianlized supplementary plot ----
plot_list <- vis_grid_clus(
  spe,
  clustervar = "discard",
  colors = c("FALSE" = "grey90", "TRUE" = "deeppink"),
  sort_clust = FALSE,
  sample_order = unique(spe$sample_id) |> sort(),
  spatial = FALSE,
  point_size = 2,
  return_plots = TRUE
)

## Plot legends ----
# pdf(
#   file = here::here(
#     "plots", "02_visium_qc",
#     paste0("vis_clus_sample_aware_low_lib_size_sfigur_legend.pdf")
#   ), height = 8, width = 8
# )
# print(cowplot::plot_grid(plotlist = plot_list[1]))
# dev.off()

## Spot Plot ----
# plot_list <- lapply(plot_list, function(p) {
#   p +
#     ggplot2::theme(
#       legend.position = "none",
#       plot.title = ggplot2::element_text(size = 40)
#     )
# })

pdf(
  file = here::here(
    "plots", "02_visium_qc",
    paste0("vis_clus_sample_aware_low_lib_size_sfigur.pdf")
  ),
  height = 6 * 8, width = 6 * 8 - 1
)
# Plot png for place holder in google drive.
# png(
#   file = here::here(
#     "plots", "02_visium_qc",
#     paste0("vis_clus_sample_aware_low_lib_size_sfigur.png")
#   ),
#   height = 5 * 8, width = 8 * 8
# )

# Plot ntc and scz in separate groups
idx_list <- spe$dx |>
  unique() |>
  set_names() |>
  imap(~ str_detect(names(plot_list), .x))

for (.idx in idx_list) {
  print(
    ggpubr::ggarrange(
      plotlist = plot_list[.idx],
      ncol = 6, nrow = 6,
      common.legend = TRUE
    )
  )
}
dev.off()

## Remove out-tissue spots ----
### Categorized outlier  (Deprecated) ----
# vis_grid_clus(
#   spe,
#   clustervar = "all_outlier",
#   spatial = FALSE,
#   pdf_file = here(
#     "plots/02_visium_qc/spot_plot_outliers.pdf"
#   ),
#   sample_order = unique(spe$sample_id) |> sort(),
#   height = 1056,
#   width = 816,
#   point_size = 0.8
# )



# Calculate per-sample statistics ----
col_dat <- colData(spe) |> data.frame()

remove_df <- col_dat |>
  group_by(sample_id) |>
  summarize(
    n_remove = sum(remove)
  )


# Total number of spots discarded
sum(remove_df$n_remove)
# [1] 3001
# Prop of spots discarded
sum(remove_df$n_remove) / nrow(col_dat)
# [1] 0.01061148

# Range of number of spots discarded in sample
remove_df$n_remove |> range()
# [1]    4 1077

# Median number of spots discarded
remove_df$n_remove |> median()
# [1] 14


# Session Info ----
session_info()
