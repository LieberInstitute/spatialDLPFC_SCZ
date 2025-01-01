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

# Plots ----


# Per-sample distribution
# sample_id_by_dx <- metadata(spe)$dx_df |>
#   arrange(dx, sample_id)



# TODO: this might be change to scatter more ----
png(
  here("plots/02_visium_qc/post_qc_metric.png"),
  width = 816, height = 1056, units = "px"
  # width = 7.5, height = 9, units = "in"
)
# gridExtra::grid.arrange(
ggpubr::ggarrange(
  # cowplot::plot_grid(
  ## sum_umi
  plotColData(spe,
    x = "sample_id", y = "sum_umi",
    color = "discard",
    point_size = 0.5
  ) +
    scale_y_log10() +
    facet_wrap(~ spe$dx, scales = "free_x") +
    # scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    labs(
      title = "Total UMIs",
      y = "library size",
      x = "sample"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      text = element_text(size = 20),
      axis.text = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.justification = "center"
    ),
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)),

  # ## sum_gene
  plotColData(spe,
    x = "sample_id", y = "sum_gene",
    color = "discard",
    point_size = 0.5
  ) +
    facet_wrap(~ spe$dx, scales = "free_x") +
    scale_y_log10() +
    # scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    labs(
      title = "Detected genes",
      y = "detected genes",
      x = "sample"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      text = element_text(size = 20),
      axis.text = element_text(size = 20)
    ),
  # # theme(axis.text.x = element_text(angle = 90, hjust = 1)),

  # ## Mito_ratio
  plotColData(spe,
    x = "sample_id", y = "expr_chrM_ratio",
    color = "discard",
    point_size = 0.5
  ) +
    facet_wrap(~ spe$dx, scales = "free_x") +
    # scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    # scale_y_log10() +
    labs(
      title = "Mitochondrial gene expression percentage",
      y = "mito gene prop.",
      x = "sample"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text = element_text(size = 20),
      text = element_text(size = 20),
      axis.text.x = element_text(
        size = 15, angle = 90, hjust = 1
      )
    ),
  ncol = 1,
  # heights = c(2, 2, 2.3),
  # widths = c(1, 1, 0.7),
  align = "v",
  legend = "bottom",
  common.legend = TRUE
)
dev.off()



## Deprecated plots ----
### Raw data spots ----
# raw_col_df <- colData(spe) |> data.frame()

# in_tissue_spots_df <- raw_col_df |>
#   group_by(sample_id) |>
#   summarize(n_in_tissue = sum(in_tissue)) |>
#   left_join(
#     metadata(spe)$dx_df,
#     by = "sample_id"
#   )

# in_tissue_spots_df |>
#   ggplot(aes(x = dx, y = n_in_tissue, color = dx)) +
#   geom_violin() +
#   geom_boxplot() #+
# #  geom_jitter()
### Outlier spots ----
# outlier_spots_df <- raw_col_df |>
#   filter(
#     in_tissue == TRUE,
#     outlier == TRUE
#   ) |>
#   group_by(sample_id) |>
#   summarize(n = n()) |>
#   left_join(
#     metadata(spe)$dx_df,
#     by = "sample_id"
#   )

# outlier_spots_df |>
#   ggplot(aes(x = dx, y = n, color = dx)) +
#   geom_violin() +
#   geom_boxplot() +
#   labs(
#     title = "# of Ouliters"
#   )
### Post-qc spots ----
# post_qc_spots_df <- raw_col_df |>
#   filter(
#     in_tissue == TRUE,
#     outlier == FALSE
#   ) |>
#   group_by(sample_id) |>
#   summarize(n = n()) |>
#   left_join(
#     metadata(spe)$dx_df,
#     by = "sample_id"
#   )

# post_qc_spots_df |>
#   ggplot(aes(x = dx, y = n, color = dx)) +
#   geom_violin() +
#   geom_boxplot() +
#   labs(
#     title = "# of Kept Spots"
#   )

# Session Info ----
sessioninfo::session_info()
