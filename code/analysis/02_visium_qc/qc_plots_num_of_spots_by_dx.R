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
# Load Data ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

## Load outlier keys ----
outlier_spots_key <- readRDS(here("processed-data/visium_qc/outlier_key_ngb_36.rds"))

## Merge data together ----
# TODO: edit this section to make it more compatible with
spe$outlier <- FALSE
spe$outlier[spe$key %in% outlier_spots_key] <- TRUE


# Remove out_tissue_spots ----
spe <- spe[, spe$in_tissue==TRUE]

# Plots ----
## Raw data spots ----
raw_col_df <- colData(spe) |> data.frame()

in_tissue_spots_df <- raw_col_df |>
  group_by(sample_id) |>
  summarize(n_in_tissue = sum(in_tissue)) |>
  left_join(
    metadata(spe)$dx_df,
    by = "sample_id"
  )

in_tissue_spots_df |>
  ggplot(aes(x = dx, y = n_in_tissue, color = dx)) +
  geom_violin() +
  geom_boxplot() #+
#  geom_jitter()


## Outlier spots ----
outlier_spots_df <- raw_col_df |>
  filter(
    in_tissue == TRUE,
    outlier == TRUE
  ) |>
  group_by(sample_id) |>
  summarize(n = n()) |>
  left_join(
    metadata(spe)$dx_df,
    by = "sample_id"
  )

outlier_spots_df |>
  ggplot(aes(x = dx, y = n, color = dx)) +
  geom_violin() +
  geom_boxplot() +
  labs(
    title = "# of Ouliters"
  )


## Post-qc spots ----
post_qc_spots_df <- raw_col_df |>
  filter(
    in_tissue == TRUE,
    outlier == FALSE
  ) |>
  group_by(sample_id) |>
  summarize(n = n()) |>
  left_join(
    metadata(spe)$dx_df,
    by = "sample_id"
  )

post_qc_spots_df |>
  ggplot(aes(x = dx, y = n, color = dx)) +
  geom_violin() +
  geom_boxplot() +
  labs(
    title = "# of Kept Spots"
  )

# Per-sample distribution
sample_id_by_dx <- metadata(spe)$dx_df |>
  arrange(dx, sample_id)



png(
  here("plots/02_visium_qc/post_qc_metric.png"),
  width = 816, height = 352, units = "px"
  # width = 8.5, height = 11/3, units = "in"
)
gridExtra::grid.arrange(
  ## sum_umi
  plotColData(spe,
    x = "sample_id", y = "sum_umi",
    color = "outlier",
    point_size = 0.5
  ) +
    scale_y_log10() +
    scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)),

  ## sum_gene
  plotColData(spe,
    x = "sample_id", y = "sum_gene",
    color = "outlier",
    point_size = 0.5
  ) +
    scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    scale_y_log10() +
    scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)),

  ## Mito_ratio
  plotColData(spe,
    x = "sample_id", y = "expr_chrM_ratio",
    color = "outlier",
    point_size = 0.5
  ) +
    scale_x_discrete(limits = sample_id_by_dx$sample_id) +
    # scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)),
  ncol = 1
)
dev.off()

# Session Info ----
sessioninfo::session_info()
