# Load library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
  library(escheR)
  library(ggpubr)
  library(cowplot)
  library(sessioninfo)
})

# Load data ----
## Load raw spe -----
# spe <- readRDS(
#   here::here(
#     "processed-data", "rds", "01_build_spe",
#     "raw_spe_wo_SPG_N63.rds"
#   )
# )
#
# ## Fetch demo info
# spe$dx <- metadata(spe)$dx_df$dx[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )
# ]
#
# spe$brnum <- metadata(spe)$dx_df$subject[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )
# ]
#
# # Create sample_label
# spe$sample_label <- paste0(
#   spe$brnum, "_", toupper(spe$dx)
# )
#
# ## Load outlier df ----
# tot_outlier_df <- readRDS(
#   here(
#     "processed-data/rds/02_visium_qc",
#     "combined_outlier_df.rds"
#   )
# )
#
# ## assemble outlier into spe object
# spe$artifact <- tot_outlier_df[spe$key, "artifact"]
# spe <- spe[, spe$in_tissue == TRUE]
#
# # Not useful
# # TODO: remove
# # spe$all_outlier <- tot_outlier_df[spe$key, "all_outlier"]
# # spe$discard <- tot_outlier_df[spe$key, "remove"]
#
# # Subset to compromised samples ----
# ## Compromised sample id ----
# vec_comp_samples <- c(
#   "V13M06-279_A1",
#   "V13M06-282_B1"
# )
#
# ## Subset to compromised samples ----
# sub_spe <- spe[, spe$sample_id %in% vec_comp_samples]
# write_rds(
#   sub_spe,
#   here(
#     "processed-data/rds/02_visium_qc",
#     "artifact_sample_spe.rds"
#   )
# )

sub_spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "artifact_sample_spe.rds"
  )
)

stopifnot(length(unique(sub_spe$sample_id)) == 2)

# Make plots ----
## UMI plots -----
umi_scale_min <- min(sub_spe$sum_umi)
umi_scale_max <- max(sub_spe$sum_umi)

p_UMI_list <- sub_spe$sample_label |>
  unique() |>
  set_names() |>
  map(.f = function(.sample) {
    # browser()
    one_spe <- sub_spe[, sub_spe$sample_label == .sample]
    make_escheR(one_spe) |>
      add_fill(
        var = "sum_umi",
        point_size = 1 # Note fix this
      ) +
      scale_fill_viridis_c(
        limits = c(umi_scale_min, umi_scale_max),
        trans = "log10",
        name = "Total UMI",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      ylab(
        .sample
      ) +
      guides(
        fill = guide_colourbar(
          # Make legend title above and centered
          title.position = "top", title.hjust = 0.5,
          # Adjust the barwidth to make the legend color scale bar longer)
          barwidth = unit(3, "cm")
        )
      ) +
      # Make legend horizontal on top
      theme(
        # legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12), # set legend font size
        legend.position = "none", # hide legend for individual plot
        # add panel border
        panel.border = element_rect(
          color = "black", fill = NA, linewidth = 1
        ),
        # set font size to 12 to be readable
        text = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90, margin = margin(r = 5)),
        plot.margin = margin(b = 10)
      )
  })


# Extract the legend from the first plot
p_UMI_legend <- ggpubr::get_legend(
  p_UMI_list[[1]] + theme(legend.position = "top"),
  position = "top"
) |> as_ggplot()

# Combine the legend and the plots
umi_plot <- plot_grid(
  p_UMI_legend,
  plot_grid(
    plotlist = p_UMI_list, ncol = 1, nrow = 2, align = "v"
  ),
  ncol = 1,
  rel_heights = c(0.12, 1), # Adjust the relative heights as needed,
  align = "v"
)
# umi_plot

## Unique genes -----
gene_scale_min <- min(sub_spe$sum_gene)
gene_scale_max <- max(sub_spe$sum_gene)

p_gene_list <- sub_spe$sample_label |>
  unique() |>
  set_names() |>
  map(.f = function(.sample) {
    one_spe <- sub_spe[, sub_spe$sample_label == .sample]
    make_escheR(one_spe) |>
      add_fill(
        var = "sum_gene",
        point_size = 1 # Note fix this
      ) +
      scale_fill_viridis_c(
        limits = c(gene_scale_min, gene_scale_max),
        trans = "log10",
        name = "Unique Genes",
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      guides(
        fill = guide_colourbar(
          # Make legend title above and centered
          title.position = "top", title.hjust = 0.5,
          # Adjust the barwidth to make the legend color scale bar longer)
          barwidth = unit(3, "cm")
        )
      ) +
      # Make legend horizontal on top
      theme(
        # legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12), # set legend font size
        legend.position = "none", # hide legend for individual plot
        # add panel border
        panel.border = element_rect(
          color = "black", fill = NA, linewidth = 1
        ),
        # set font size to 12 to be readable
        text = element_text(size = 12),
        # axis.title.y = element_text(size = 12, angle = 90, margin = margin(r = 10)),
        plot.margin = margin(b = 10, l = 17)
      )
  })


# Extract the legend from the first plot
p_gene_legend <- ggpubr::get_legend(
  p_gene_list[[1]] + theme(legend.position = "top"),
  position = "top"
) |> as_ggplot()

# Combine the legend and the plots
gene_plot <- plot_grid(
  p_gene_legend,
  plot_grid(
    plotlist = p_gene_list, ncol = 1, nrow = 2, align = "v"
  ),
  ncol = 1,
  rel_heights = c(0.12, 1), # Adjust the relative heights as needed,
  align = "v"
)
# gene_plot

### Mitocondria Ratio -----
mito_scale_min <- min(sub_spe$expr_chrM_ratio)
mito_scale_max <- max(sub_spe$expr_chrM_ratio)

p_mito_list <- sub_spe$sample_label |>
  unique() |>
  set_names() |>
  map(.f = function(.sample) {
    one_spe <- sub_spe[, sub_spe$sample_label == .sample]
    make_escheR(one_spe) |>
      add_fill(
        var = "expr_chrM_ratio",
        point_size = 1 # Note fix this
      ) +
      scale_fill_viridis_c(
        limits = c(mito_scale_min, mito_scale_max),
        name = "Mito. Raio"
      ) +
      guides(
        fill = guide_colourbar(
          # Make legend title above and centered
          title.position = "top", title.hjust = 0.5,
          # Adjust the barwidth to make the legend color scale bar longer)
          barwidth = unit(3, "cm")
        )
      ) +
      # Make legend horizontal on top
      theme(
        # legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12), # set legend font size
        legend.position = "none", # hide legend for individual plot
        # add panel border
        panel.border = element_rect(
          color = "black", fill = NA, linewidth = 1
        ),
        # set font size to 12 to be readable
        text = element_text(size = 12),
        # axis.title.y = element_text(size = 12, angle = 90, margin = margin(r = 10)),
        plot.margin = margin(b = 10, l = 17)
      )
  })


# Extract the legend from the first plot
p_mito_legend <- ggpubr::get_legend(
  p_mito_list[[1]] + theme(legend.position = "top"),
  position = "top"
) |> as_ggplot()

# Combine the legend and the plots
mito_plot <- plot_grid(
  p_mito_legend,
  plot_grid(
    plotlist = p_mito_list, ncol = 1, nrow = 2, align = "v"
  ),
  ncol = 1,
  rel_heights = c(0.12, 1), # Adjust the relative heights as needed,
  align = "v"
)
# mito_plot


## Artifacts plot ----
p_artifact_list <- sub_spe$sample_label |>
  unique() |>
  set_names() |>
  map(.f = function(.sample) {
    one_spe <- sub_spe[, sub_spe$sample_label == .sample]
    ret_plot <- make_escheR(one_spe) |>
      add_fill(var = "artifact", point_size = 1)

    # organize legends
    ret_plot <- ret_plot +
      scale_fill_manual(
        name = "Regional Artifect",
        limits = c("TRUE", "FALSE"), # Show discard first
        values = c(
          "TRUE" = "deeppink",
          "FALSE" = "grey90"
        ),
        labels = c(
          "TRUE" = "Yes",
          "FALSE" = "No"
        )
      ) +
      # Make legend circle more obvious
      guides(
        fill = guide_legend(
          title.position = "top", title.hjust = 0.5,
          override.aes = list(size = 4)
        )
      )

    # organize appearance
    ret_plot <- ret_plot +
      theme(
        # add panel border
        panel.border = element_rect(
          color = "black", fill = NA, linewidth = 1
        ),
        # set font size to 12 to be readable
        text = element_text(size = 12),
        legend.title = element_text(size = 12), # set legend font size
        legend.position = "none", # hide legend for individual plot
        plot.margin = margin(b = 10, l = 17)
      )

    ret_plot
  })

p_artifact_legend <- ggpubr::get_legend(
  p_artifact_list[[1]] + theme(legend.position = "top"),
  position = "top"
) |> as_ggplot()

# Combine the legend and the plots
artifact_plot <- plot_grid(
  p_artifact_legend,
  plot_grid(
    plotlist = p_artifact_list, ncol = 1, nrow = 2, align = "v"
  ),
  ncol = 1,
  rel_heights = c(0.12, 1), # Adjust the relative heights as needed,
  align = "v"
)

# Create Panels ----
## Save plot as pds ----
pdf(
  here(
    "plots/02_visium_qc",
    "qc_regional_artifact.pdf"
  ),
  height = 5,
  width = 8
)
plot_grid(
  # Legend row
  plot_grid(
    p_UMI_legend, p_gene_legend, p_mito_legend, p_artifact_legend,
    nrow = 1, ncol = 4
  ),
  # two rows of spot plots
  plot_grid(
    plotlist = c(
      p_UMI_list,
      p_gene_list,
      p_mito_list,
      p_artifact_list
    ),
    nrow = 2,
    ncol = 4,
    byrow = FALSE,
    align = "h",
    axis = "b"
  ),
  nrow = 2, ncol = 1,
  rel_heights = c(0.2, 1)
) |> print()
dev.off()

# Session info ----
session_info()
# ─ Session info ─────────────────────────────────────────────
# setting  value
# version  R version 4.4.1 (2024-06-14)
# os       macOS Sonoma 14.6.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2025-01-07
# rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
# pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# abind                  1.4-8    2024-09-12 [1] CRAN (R 4.4.1)
# backports              1.5.0    2024-05-23 [1] CRAN (R 4.4.0)
# Biobase              * 2.64.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# broom                  1.0.7    2024-09-26 [1] CRAN (R 4.4.1)
# car                    3.1-3    2024-09-27 [1] CRAN (R 4.4.1)
# carData                3.0-5    2022-01-06 [1] CRAN (R 4.4.0)
# cli                    3.6.3    2024-06-21 [1] CRAN (R 4.4.0)
# colorspace             2.1-1    2024-07-26 [1] CRAN (R 4.4.0)
# cowplot              * 1.1.3    2024-01-22 [1] CRAN (R 4.4.0)
# crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.4.0)
# DelayedArray           0.30.1   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
# dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.4.0)
# escheR               * 1.4.0    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
# fansi                  1.0.6    2023-12-08 [1] CRAN (R 4.4.0)
# farver                 2.1.2    2024-05-13 [1] CRAN (R 4.4.0)
# forcats              * 1.0.0    2023-01-29 [1] CRAN (R 4.4.0)
# Formula                1.2-5    2023-02-24 [1] CRAN (R 4.4.0)
# generics               0.1.3    2022-07-05 [1] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.40.1   2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12   2024-08-05 [1] Bioconductor
# GenomicRanges        * 1.56.2   2024-10-09 [1] Bioconductor 3.19 (R 4.4.1)
# ggplot2              * 3.5.1    2024-04-23 [1] CRAN (R 4.4.0)
# ggpubr               * 0.6.0    2023-02-10 [1] CRAN (R 4.4.0)
# ggsignif               0.6.4    2022-10-13 [1] CRAN (R 4.4.0)
# glue                   1.8.0    2024-09-30 [1] CRAN (R 4.4.1)
# gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.4.1)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.4.0)
# hms                    1.1.3    2023-03-21 [1] CRAN (R 4.4.0)
# httr                   1.4.7    2023-08-15 [1] CRAN (R 4.4.0)
# IRanges              * 2.38.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
# jsonlite               1.8.9    2024-09-20 [1] CRAN (R 4.4.1)
# labeling               0.4.3    2023-08-29 [1] CRAN (R 4.4.0)
# lattice                0.22-6   2024-03-20 [1] CRAN (R 4.4.1)
# lifecycle              1.0.4    2023-11-07 [1] CRAN (R 4.4.0)
# lubridate            * 1.9.3    2023-09-27 [1] CRAN (R 4.4.0)
# magick                 2.8.5    2024-09-20 [1] CRAN (R 4.4.1)
# magrittr               2.0.3    2022-03-30 [1] CRAN (R 4.4.0)
# Matrix                 1.7-1    2024-10-18 [1] CRAN (R 4.4.1)
# MatrixGenerics       * 1.16.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.4.1    2024-09-08 [1] CRAN (R 4.4.1)
# munsell                0.5.1    2024-04-01 [1] CRAN (R 4.4.0)
# pillar                 1.9.0    2023-03-22 [1] CRAN (R 4.4.0)
# pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.4.0)
# purrr                * 1.0.2    2023-08-10 [1] CRAN (R 4.4.0)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.4.0)
# Rcpp                   1.0.13-1 2024-11-02 [1] CRAN (R 4.4.1)
# readr                * 2.1.5    2024-01-10 [1] CRAN (R 4.4.0)
# rjson                  0.2.23   2024-09-16 [1] CRAN (R 4.4.1)
# rlang                  1.1.4    2024-06-04 [1] CRAN (R 4.4.0)
# rprojroot              2.0.4    2023-11-05 [1] CRAN (R 4.4.0)
# rstatix                0.7.2    2023-02-01 [1] CRAN (R 4.4.0)
# rstudioapi             0.17.1   2024-10-22 [1] CRAN (R 4.4.1)
# S4Arrays               1.4.1    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
# scales                 1.3.0    2023-11-28 [1] CRAN (R 4.4.0)
# sessioninfo          * 1.2.2    2021-12-06 [1] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# SparseArray            1.4.8    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
# SpatialExperiment    * 1.14.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# stringi                1.8.4    2024-05-06 [1] CRAN (R 4.4.0)
# stringr              * 1.5.1    2023-11-14 [1] CRAN (R 4.4.0)
# SummarizedExperiment * 1.34.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# tibble               * 3.2.1    2023-03-20 [1] CRAN (R 4.4.0)
# tidyr                * 1.3.1    2024-01-24 [1] CRAN (R 4.4.0)
# tidyselect             1.2.1    2024-03-11 [1] CRAN (R 4.4.0)
# tidyverse            * 2.0.0    2023-02-22 [1] CRAN (R 4.4.0)
# timechange             0.3.0    2024-01-18 [1] CRAN (R 4.4.0)
# tzdb                   0.4.0    2023-05-12 [1] CRAN (R 4.4.0)
# UCSC.utils             1.0.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
# utf8                   1.2.4    2023-10-22 [1] CRAN (R 4.4.0)
# vctrs                  0.6.5    2023-12-01 [1] CRAN (R 4.4.0)
# viridisLite            0.4.2    2023-05-02 [1] CRAN (R 4.4.0)
# withr                  3.0.2    2024-10-28 [1] CRAN (R 4.4.1)
# XVector                0.44.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# zlibbioc               1.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
