# Load library -----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(cowplot)
  library(sessioninfo)

  # library(scater)
})

# Load Data -----
## Raw spe -----
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

spe$DX <- toupper(spe$dx)

spe$brnum <- metadata(spe)$DX_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$DX_df$sample_id
  )
]

spe$sample_label <- paste0(
  spe$brnum, "_", spe$DX
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


# Creat plots -----
## Set seed ----
# note: geom_jitter use a randomized algorithm
set.seed(20240110)

## Discard spot plots -----
### Create discard data df -----
discard_col_df <- colData(spe[, spe$discard == TRUE]) |>
  as.data.frame()

### Prop of spot dicarded ----
sample_lvl_summary_df <- colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(
    sample_label = unique(sample_label),
    DX = unique(DX),
    n_total = n(),
    n_discard = sum(discard),
    prop_discard = n_discard / n_total,
  )

discard_prop_p <- sample_lvl_summary_df |>
  ggplot(aes(x = DX, y = prop_discard)) +
  geom_boxplot(aes(color = DX)) +
  geom_jitter(size = 0.5, alpha = 0.3) +
  theme_light() +
  labs(
    title = "Discarded Spots",
    y = "Prop. of Spots"
  ) +
  scale_color_manual(
    values = c(
      "NTC" = "blue",
      "SCZ" = "red"
    )
  )

discard_prop_p
# NOTE: Optimally, Boyi would want to label the two samples that
#       contains regional artifacts in the plot.

### UMI ----
discard_umi_p <- discard_col_df |>
  ggplot() +
  stat_ecdf(
    aes(x = sum_umi, group = sample_id, color = DX),
    alpha = 0.5,
    geom = "step"
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_minimal() +
  labs(
    title = "Total UMI",
    x = "Total UMI",
    y = "Cumulative Probability",
    color = "Diagnosis"
  )
discard_umi_p

### Unique genes ----
discard_gene_p <- discard_col_df |>
  ggplot() +
  stat_ecdf(
    aes(x = sum_gene, group = sample_id, color = DX),
    alpha = 0.5,
    geom = "step"
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_minimal() +
  labs(
    title = "Unique Genes",
    x = "Unique Genes",
    y = "Cumulative Probability",
    color = "Diagnosis"
  )
discard_gene_p

### Mito Ratio ----
discard_mito_p <- discard_col_df |>
  ggplot() +
  stat_ecdf(
    aes(x = expr_chrM_ratio, group = sample_id, color = DX),
    alpha = 0.5,
    geom = "step"
  ) +
  theme_minimal() +
  labs(
    title = "Mit. Ratio",
    x = "Mit. Ratio",
    y = "Cumulative Probability",
    color = "Diagnosis"
  )
discard_mito_p

## Kept spot plots -----
### Create kept data df -----
kept_col_df <- colData(spe[, spe$discard != TRUE]) |>
  data.frame()

### Numb of kept spots ----
kept_n_plot <- kept_col_df |>
  group_by(sample_id) |>
  summarize(
    sample_label = unique(sample_label),
    DX = unique(DX),
    n = n()
  ) |>
  ggplot(aes(x = DX, y = n)) +
  geom_boxplot(aes(color = DX)) +
  geom_jitter(size = 0.5, alpha = 0.3) +
  labs(
    title = "Kept Spots",
    y = "Number of Spots"
  ) +
  theme_light()
kept_n_plot

### UMI ----
kept_umi_p <- kept_col_df |>
  ggplot() +
  stat_ecdf(
    aes(x = sum_umi, group = sample_id, color = DX),
    alpha = 0.5,
    geom = "step"
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_minimal() +
  labs(
    title = "Total UMI",
    x = "UMI Counts",
    y = "Cumulative Probability",
    color = "Diagnosis"
  )
kept_umi_p

### Unique genes ----
kept_gene_p <- kept_col_df |>
  ggplot() +
  stat_ecdf(
    aes(x = sum_gene, group = sample_id, color = DX),
    alpha = 0.5,
    geom = "step"
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_minimal() +
  labs(
    title = "Unique Genes",
    x = "Unique Genes",
    y = "Cumulative Probability",
    color = "Diagnosis"
  )
kept_gene_p

## Mito Ratio ----
kept_mito_p <- kept_col_df |>
  ggplot() +
  stat_ecdf(
    aes(x = expr_chrM_ratio, group = sample_id, color = DX),
    alpha = 0.5,
    geom = "step"
  ) +
  theme_minimal() +
  labs(
    title = "Mito. Ratio",
    x = "Mito. Ratio",
    y = "Cumulative Probability",
    color = "Diagnosis"
  )
kept_gene_p

# Create paneled plot -----
## Create plot list ----
raw_p_list <- list(
  # discard spot plots
  discard_prop_p, discard_umi_p, discard_gene_p, discard_mito_p,
  # kept spot plots
  kept_n_plot, kept_umi_p, kept_gene_p, kept_mito_p
)

## Adjust theme and labels ----
ret_p_list <- raw_p_list |> map(.f = function(.p) {
  .p +
    scale_color_manual(
      values = c(
        "NTC" = "blue",
        "SCZ" = "red"
      )
    ) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(t = 5, b = 5, l = 5, r = 20)
    )
})

## Create Planel ----

panel_p <- cowplot::plot_grid(
  plotlist = ret_p_list,
  ncol = 2, nrow = 4,
  byrow = FALSE
)
panel_p

## Save plot ----
ggsave(
  here(
    "plots/02_visium_qc",
    "qc_plot_CDF_discard_kept_spots.pdf"
  ),
  plot = panel_p,
  width = 7.5, height = 10, units = "in"
)

# Session Info ----
session_info()
# ─ Session info ─────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.1 (2024-06-14)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-01-10
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ─────────────────────────────────────────────────────────────
#  package              * version  date (UTC) lib source
#  abind                  1.4-8    2024-09-12 [1] CRAN (R 4.4.1)
#  beachmat               2.20.0   2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  beeswarm               0.4.0    2021-06-01 [1] CRAN (R 4.4.0)
#  Biobase              * 2.64.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocGenerics         * 0.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocNeighbors          1.22.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocParallel           1.38.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocSingular           1.20.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  cli                    3.6.3    2024-06-21 [1] CRAN (R 4.4.0)
#  codetools              0.2-20   2024-03-31 [1] CRAN (R 4.4.1)
#  colorspace             2.1-1    2024-07-26 [1] CRAN (R 4.4.0)
#  cowplot              * 1.1.3    2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.4.0)
#  DelayedArray           0.30.1   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  DelayedMatrixStats     1.26.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.4.0)
#  fansi                  1.0.6    2023-12-08 [1] CRAN (R 4.4.0)
#  farver                 2.1.2    2024-05-13 [1] CRAN (R 4.4.0)
#  forcats              * 1.0.0    2023-01-29 [1] CRAN (R 4.4.0)
#  generics               0.1.3    2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.40.1   2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomeInfoDbData       1.2.12   2024-08-05 [1] Bioconductor
#  GenomicRanges        * 1.56.2   2024-10-09 [1] Bioconductor 3.19 (R 4.4.1)
#  ggbeeswarm             0.7.2    2023-04-29 [1] CRAN (R 4.4.0)
#  ggplot2              * 3.5.1    2024-04-23 [1] CRAN (R 4.4.0)
#  ggrepel                0.9.6    2024-09-07 [1] CRAN (R 4.4.1)
#  glue                   1.8.0    2024-09-30 [1] CRAN (R 4.4.1)
#  gridExtra              2.3      2017-09-09 [1] CRAN (R 4.4.0)
#  gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.4.1)
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3    2023-03-21 [1] CRAN (R 4.4.0)
#  httpgd                 2.0.2    2024-06-05 [1] CRAN (R 4.4.0)
#  httr                   1.4.7    2023-08-15 [1] CRAN (R 4.4.0)
#  IRanges              * 2.38.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  irlba                  2.3.5.1  2022-10-03 [1] CRAN (R 4.4.0)
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
#  ragg                   1.3.3    2024-09-11 [1] CRAN (R 4.4.1)
#  Rcpp                   1.0.13-1 2024-11-02 [1] CRAN (R 4.4.1)
#  readr                * 2.1.5    2024-01-10 [1] CRAN (R 4.4.0)
#  rjson                  0.2.23   2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.4    2024-06-04 [1] CRAN (R 4.4.0)
#  rprojroot              2.0.4    2023-11-05 [1] CRAN (R 4.4.0)
#  rsvd                   1.0.5    2021-04-16 [1] CRAN (R 4.4.0)
#  S4Arrays               1.4.1    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Vectors            * 0.42.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  ScaledMatrix           1.12.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  scales                 1.3.0    2023-11-28 [1] CRAN (R 4.4.0)
#  scater               * 1.32.1   2024-07-21 [1] Bioconductor 3.19 (R 4.4.1)
#  scuttle              * 1.14.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sessioninfo          * 1.2.2    2021-12-06 [1] CRAN (R 4.4.0)
#  SingleCellExperiment * 1.26.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SparseArray            1.4.8    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sparseMatrixStats      1.16.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SpatialExperiment    * 1.14.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  stringi                1.8.4    2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1    2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.34.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  systemfonts            1.1.0    2024-05-15 [1] CRAN (R 4.4.0)
#  textshaping            0.4.0    2024-05-24 [1] CRAN (R 4.4.0)
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
#  vipor                  0.4.7    2023-12-18 [1] CRAN (R 4.4.0)
#  viridis                0.6.5    2024-01-29 [1] CRAN (R 4.4.0)
#  viridisLite            0.4.2    2023-05-02 [1] CRAN (R 4.4.0)
#  withr                  3.0.2    2024-10-28 [1] CRAN (R 4.4.1)
#  XVector                0.44.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  zlibbioc               1.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

# ────────────────────────────────────────────────────────────────────────