# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(readxl)
  library(escheR)
  library(ggpubr)
  library(here)
  library(sessioninfo)
})

# Load Data ----
## Load Spe ----
raw_spe <- readRDS(
  here(
    # "processed-data/rds/spatial_cluster",
    # "PRECAST",
    # "spe_wo_spg_N63_PRECAST.rds"
    "processed-data/rds/01_build_spe/fnl_spe_kept_spots_only.rds"
  )
)

# error prevention
stopifnot(all(raw_spe$in_tissue))

# create sample_label
raw_spe$sample_label <- paste0(
  raw_spe$brnum, "_", toupper(raw_spe$dx)
)

# ## Load PRECAST df ----
# PRECAST_df <- readRDS(
#   here(
#     "processed-data/rds/spatial_cluster",
#     "PRECAST",
#     "test_clus_label_df_semi_inform_k_2-16.rds"
#   )
# )

## Merge PRECAST df ----
# precast_vars <- grep(
#   "^PRECAST_", colnames(PRECAST_df),
#   value = TRUE
# )

spe <- raw_spe#[, raw_spe$key %in% PRECAST_df$key]
# raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
# col_data_df <- PRECAST_df |>
#   right_join(
#     colData(spe) |> data.frame(),
#     by = c("key"),
#     relationship = "one-to-one"
#   )
# rownames(col_data_df) <- colnames(spe)
# colData(spe) <- DataFrame(col_data_df)

spe$PRECAST_07 <- spe$spd_label |> as.character() |> str_sub(start = 1, end = 5) |> factor()
# error prevention
# stopifnot(is.character(spe$PRECAST_07))

# Subset to 24 Xenium Samples ----
## Find Xenium Samples ----
xn_samples_name <- read_xlsx(
  here::here("code/xenium_panel_design/Xenium_DonorList_Edit.xlsx"),
  col_names = FALSE
)[, 1:2] |> unlist()
xn_samples_name <- xn_samples_name[!is.na(xn_samples_name)]
stopifnot(length(xn_samples_name) == 24)

## SPE with only Xenium Samples ----
spe <- spe[, spe$brnum %in% xn_samples_name]

# SPE info
# r$> spe
# class: SpatialExperiment
# dim: 36601 109353

# error prevention
stopifnot(length(spe$sample_id |> unique()) == 24)
stopifnot(ncol(spe) != 0)


# Create spot plot with escheR ----

## Create a list of escheR plots for samples ----
spot_plot_list <- spe$sample_label |>
  unique() |>
  set_names() |>
  map(.f = function(.smp) {
    sub_spe <- spe[, spe$sample_label == .smp]

    ret_plot <- sub_spe |>
      make_escheR() |>
      add_fill(var = "PRECAST_07", point_size = 0.8)

    # organize legends
    ret_plot <- ret_plot +
      scale_fill_manual(
        name = "Spatial Domain",
        values = set_names(
          Polychrome::palette36.colors(7)[seq.int(7)],
          unique(spe$PRECAST_07) |> sort()
        )
      ) +
      # Make legend circle more obvious
      guides(fill = guide_legend(override.aes = list(size = 5)))

    # organize appearance
    ret_plot <- ret_plot +
      labs(title = .smp) +
      theme(
        # adjust title
        plot.title = element_text(
          hjust = 0.5, # center title
          vjust = 2, # raise title position
          size = 12 # set font size
        ),
        # add panel border
        panel.border = element_rect(
          color = "black", fill = NA, linewidth = 1
        ),
        # set font size to 12 to be readable
        text = element_text(size = 12),
        legend.text = element_text(size = 12), # set legend font size
        legend.position = "none" # hide legend for individual plot
      )

    ret_plot
  })

# error prevention
stopifnot(
  length(names(spot_plot_list)) == length(xn_samples_name)
)

## Order samples ----
# Order samples by Dx, followed by slide_numbers
smp_order <- colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(dx, sample_id, sample_label) |>
  arrange(dx, sample_id) |>
  pull(sample_label)

## Save pdf output ----
# portrait 6 (3+3) rows and 4 columns
pdf(
  here(
    "plots/03_visium_spatial_clustering",
    "spot_plot_PRECAST07_24_samples.pdf"
  ),
  # NOTE: You need to adjust escheR point size
  # `add_fill(var = "discard", size = 1.2)`
  # if the size of the figure is changed.
  height = 11, width = 7
)
ggpubr::ggarrange(
  plotlist = spot_plot_list[smp_order],
  nrow = 6, ncol = 4,
  common.legend = TRUE,
  legend = "top"
)
dev.off()

# Session Info ----
session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.1 (2024-06-14)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-01-06
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ──────────────────────────────────────────────────────────────────────────────
#  package              * version  date (UTC) lib source
#  abind                  1.4-8    2024-09-12 [1] CRAN (R 4.4.1)
#  backports              1.5.0    2024-05-23 [1] CRAN (R 4.4.0)
#  Biobase              * 2.64.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocGenerics         * 0.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  broom                  1.0.7    2024-09-26 [1] CRAN (R 4.4.1)
#  car                    3.1-3    2024-09-27 [1] CRAN (R 4.4.1)
#  carData                3.0-5    2022-01-06 [1] CRAN (R 4.4.0)
#  cellranger             1.1.0    2016-07-27 [1] CRAN (R 4.4.0)
#  cli                    3.6.3    2024-06-21 [1] CRAN (R 4.4.0)
#  colorspace             2.1-1    2024-07-26 [1] CRAN (R 4.4.0)
#  cowplot                1.1.3    2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.4.0)
#  DelayedArray           0.30.1   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.4.0)
#  escheR               * 1.4.0    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  fansi                  1.0.6    2023-12-08 [1] CRAN (R 4.4.0)
#  farver                 2.1.2    2024-05-13 [1] CRAN (R 4.4.0)
#  forcats              * 1.0.0    2023-01-29 [1] CRAN (R 4.4.0)
#  Formula                1.2-5    2023-02-24 [1] CRAN (R 4.4.0)
#  generics               0.1.3    2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.40.1   2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomeInfoDbData       1.2.12   2024-08-05 [1] Bioconductor
#  GenomicRanges        * 1.56.2   2024-10-09 [1] Bioconductor 3.19 (R 4.4.1)
#  ggplot2              * 3.5.1    2024-04-23 [1] CRAN (R 4.4.0)
#  ggpubr                 0.6.0    2023-02-10 [1] CRAN (R 4.4.0)
#  ggsignif               0.6.4    2022-10-13 [1] CRAN (R 4.4.0)
#  glue                   1.8.0    2024-09-30 [1] CRAN (R 4.4.1)
#  gridExtra              2.3      2017-09-09 [1] CRAN (R 4.4.0)
#  gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.4.1)
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3    2023-03-21 [1] CRAN (R 4.4.0)
#  httpgd                 2.0.2    2024-06-05 [1] CRAN (R 4.4.0)
#  httr                   1.4.7    2023-08-15 [1] CRAN (R 4.4.0)
#  IRanges              * 2.38.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  jsonlite               1.8.9    2024-09-20 [1] CRAN (R 4.4.1)
#  labeling               0.4.3    2023-08-29 [1] CRAN (R 4.4.0)
#  lattice                0.22-6   2024-03-20 [1] CRAN (R 4.4.1)
#  lifecycle              1.0.4    2023-11-07 [1] CRAN (R 4.4.0)
#  lubridate            * 1.9.3    2023-09-27 [1] CRAN (R 4.4.0)
#  magick                 2.8.5    2024-09-20 [1] CRAN (R 4.4.1)
#  magrittr               2.0.3    2022-03-30 [1] CRAN (R 4.4.0)
#  Matrix                 1.7-1    2024-10-18 [1] CRAN (R 4.4.1)
#  MatrixGenerics       * 1.16.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  matrixStats          * 1.4.1    2024-09-08 [1] CRAN (R 4.4.1)
#  munsell                0.5.1    2024-04-01 [1] CRAN (R 4.4.0)
#  pillar                 1.9.0    2023-03-22 [1] CRAN (R 4.4.0)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.4.0)
#  purrr                * 1.0.2    2023-08-10 [1] CRAN (R 4.4.0)
#  R6                     2.5.1    2021-08-19 [1] CRAN (R 4.4.0)
#  Rcpp                   1.0.13-1 2024-11-02 [1] CRAN (R 4.4.1)
#  readr                * 2.1.5    2024-01-10 [1] CRAN (R 4.4.0)
#  readxl               * 1.4.3    2023-07-06 [1] CRAN (R 4.4.0)
#  rjson                  0.2.23   2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.4    2024-06-04 [1] CRAN (R 4.4.0)
#  rprojroot              2.0.4    2023-11-05 [1] CRAN (R 4.4.0)
#  rstatix                0.7.2    2023-02-01 [1] CRAN (R 4.4.0)
#  S4Arrays               1.4.1    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Vectors            * 0.42.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  scales                 1.3.0    2023-11-28 [1] CRAN (R 4.4.0)
#  sessioninfo          * 1.2.2    2021-12-06 [1] CRAN (R 4.4.0)
#  SingleCellExperiment * 1.26.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SparseArray            1.4.8    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SpatialExperiment    * 1.14.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  stringi                1.8.4    2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1    2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.34.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  systemfonts            1.1.0    2024-05-15 [1] CRAN (R 4.4.0)
#  tibble               * 3.2.1    2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                * 1.3.1    2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1    2024-03-11 [1] CRAN (R 4.4.0)
#  tidyverse            * 2.0.0    2023-02-22 [1] CRAN (R 4.4.0)
#  timechange             0.3.0    2024-01-18 [1] CRAN (R 4.4.0)
#  tzdb                   0.4.0    2023-05-12 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.0.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  unigd                  0.1.2    2024-06-05 [1] CRAN (R 4.4.0)
#  utf8                   1.2.4    2023-10-22 [1] CRAN (R 4.4.0)
#  vctrs                  0.6.5    2023-12-01 [1] CRAN (R 4.4.0)
#  withr                  3.0.2    2024-10-28 [1] CRAN (R 4.4.1)
#  XVector                0.44.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  zlibbioc               1.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library