# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(escheR)
  library(sessioninfo)
})


# Load Data ----
# note(Boyi): It would be much eaiser if Boyi created a spe object only contains
# the two representative samples and relevant colData.
raw_spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# error prevention
stopifnot(all(raw_spe$in_tissue))

# create sample_label
raw_spe$sample_label <- paste0(
  raw_spe$brnum, "_", toupper(raw_spe$dx)
)

# ## Load PRECAST df ----
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

spe <- raw_spe[, raw_spe$key %in% PRECAST_df$key]
col_data_df <- PRECAST_df |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )
rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# error prevention
stopifnot(is.character(spe$PRECAST_07))

# Subset to representative samples ----
rep_sample_id <- c("V13M06-342_D1", "V13M06-343_D1")

# error prevention
stopifnot(
  all(
    rep_sample_id %in%
      unique(spe$sample_id)
  )
)

spe <- spe[, spe$sample_id %in% rep_sample_id]

# error prevention
stopifnot(
  length(
    unique(spe$sample_id)
  ) == 2
)

# Create multi-gene spot plot ----
# genes include MBP, PCP4, SNAP25
plot_list <- unique(spe$sample_label) |>
  set_names() |>
  map(.f = function(.smp) {
    sub_spe <- spe[, spe$sample_label == .smp]

    make_escheR(sub_spe) |>
      add_fill(
        "PRECAST_07",
        point_size = 2.1
      ) +
      labs(title = .smp) +
      theme(
        # legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5),
        # panel.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )
  })


## Adjust color palette for PRECAST_07 ----
plot_list <- plot_list |> lapply(FUN = function(.p) {
  .p +
    scale_fill_manual(
      name = "Spatial Domain",
      values = set_names(
        Polychrome::palette36.colors(7)[seq.int(7)],
        unique(spe$PRECAST_07) |> sort()
      ),
      guide = guide_legend(override.aes = list(size = 7))
    ) +
    theme(
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15)
    )
})

# Combine plots with legend ----
combined_plot <- ggpubr::ggarrange(
  plotlist = plot_list,
  nrow = 1,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom" # ,
  # legend.grob = ggpubr::get_legend(legend_plot)
)

# Save combined plot ----
ggsave(
  filename = here(
    "plots/03_visium_spatial_clustering",
    "spd_spot_plot_PRECAST-07_rep_sample.pdf"
  ),
  plot = combined_plot,
  height = 5.5, width = 9.5,
  unit = "in"
)

# Session Info
session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.1 (2024-06-14)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-01-15
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ──────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-8     2024-09-12 [1] CRAN (R 4.4.1)
#  AnnotationDbi          1.66.0    2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  AnnotationHub          3.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  attempt                0.3.1     2020-05-03 [1] CRAN (R 4.4.0)
#  backports              1.5.0     2024-05-23 [1] CRAN (R 4.4.0)
#  beachmat               2.20.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.4.0)
#  benchmarkme            1.0.8     2022-06-12 [1] CRAN (R 4.4.0)
#  benchmarkmeData        1.0.4     2020-04-23 [1] CRAN (R 4.4.0)
#  Biobase              * 2.64.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocFileCache          2.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocGenerics         * 0.50.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocIO                 1.14.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocManager            1.30.25   2024-08-28 [1] CRAN (R 4.4.1)
#  BiocNeighbors          1.22.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocParallel           1.38.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocSingular           1.20.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocVersion            3.19.1    2024-04-22 [1] Bioconductor 3.19 (R 4.4.0)
#  Biostrings             2.72.1    2024-06-02 [1] Bioconductor 3.19 (R 4.4.0)
#  bit                    4.5.0.1   2024-12-03 [1] CRAN (R 4.4.1)
#  bit64                  4.5.2     2024-09-22 [1] CRAN (R 4.4.1)
#  bitops                 1.0-9     2024-10-03 [1] CRAN (R 4.4.1)
#  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.4.0)
#  broom                  1.0.7     2024-09-26 [1] CRAN (R 4.4.1)
#  bslib                  0.8.0     2024-07-29 [1] CRAN (R 4.4.0)
#  cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.4.0)
#  car                    3.1-3     2024-09-27 [1] CRAN (R 4.4.1)
#  carData                3.0-5     2022-01-06 [1] CRAN (R 4.4.0)
#  cli                    3.6.3     2024-06-21 [1] CRAN (R 4.4.0)
#  codetools              0.2-20    2024-03-31 [1] CRAN (R 4.4.1)
#  colorspace             2.1-1     2024-07-26 [1] CRAN (R 4.4.0)
#  config                 0.3.2     2023-08-30 [1] CRAN (R 4.4.0)
#  cowplot                1.1.3     2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3     2024-06-20 [1] CRAN (R 4.4.0)
#  curl                   6.1.0     2025-01-06 [1] CRAN (R 4.4.1)
#  data.table             1.16.4    2024-12-06 [1] CRAN (R 4.4.1)
#  DBI                    1.2.3     2024-06-02 [1] CRAN (R 4.4.0)
#  dbplyr                 2.5.0     2024-03-19 [1] CRAN (R 4.4.0)
#  DelayedArray           0.30.1    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  DelayedMatrixStats     1.26.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  digest                 0.6.37    2024-08-19 [1] CRAN (R 4.4.1)
#  doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.4.0)
#  dotCall64              1.2       2024-10-04 [1] CRAN (R 4.4.1)
#  dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.4.0)
#  DT                     0.33      2024-04-04 [1] CRAN (R 4.4.0)
#  edgeR                  4.2.2     2024-10-13 [1] Bioconductor 3.19 (R 4.4.1)
#  escheR               * 1.4.0     2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  ExperimentHub          2.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  farver                 2.1.2     2024-05-13 [1] CRAN (R 4.4.0)
#  fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.4.0)
#  fields                 16.3      2024-09-30 [1] CRAN (R 4.4.1)
#  filelock               1.0.3     2023-12-11 [1] CRAN (R 4.4.0)
#  forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.4.0)
#  foreach                1.5.2     2022-02-02 [1] CRAN (R 4.4.0)
#  Formula                1.2-5     2023-02-24 [1] CRAN (R 4.4.0)
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.40.1    2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomeInfoDbData       1.2.12    2024-08-05 [1] Bioconductor
#  GenomicAlignments      1.40.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomicRanges        * 1.56.2    2024-10-09 [1] Bioconductor 3.19 (R 4.4.1)
#  ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.4.0)
#  ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.4.0)
#  ggpubr               * 0.6.0     2023-02-10 [1] CRAN (R 4.4.0)
#  ggrepel                0.9.6     2024-09-07 [1] CRAN (R 4.4.1)
#  ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.4.0)
#  glue                   1.8.0     2024-09-30 [1] CRAN (R 4.4.1)
#  golem                  0.5.1     2024-08-27 [1] CRAN (R 4.4.1)
#  gridExtra              2.3       2017-09-09 [1] CRAN (R 4.4.0)
#  gtable                 0.3.6     2024-10-25 [1] CRAN (R 4.4.1)
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.4.0)
#  htmltools              0.5.8.1   2024-04-04 [1] CRAN (R 4.4.0)
#  htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.4.0)
#  httpgd                 2.0.2     2024-06-05 [1] CRAN (R 4.4.0)
#  httpuv                 1.6.15    2024-03-26 [1] CRAN (R 4.4.0)
#  httr                   1.4.7     2023-08-15 [1] CRAN (R 4.4.0)
#  IRanges              * 2.38.1    2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.4.0)
#  iterators              1.0.14    2022-02-05 [1] CRAN (R 4.4.0)
#  jquerylib              0.1.4     2021-04-26 [1] CRAN (R 4.4.0)
#  jsonlite               1.8.9     2024-09-20 [1] CRAN (R 4.4.1)
#  KEGGREST               1.44.1    2024-06-19 [1] Bioconductor 3.19 (R 4.4.0)
#  labeling               0.4.3     2023-08-29 [1] CRAN (R 4.4.0)
#  later                  1.4.1     2024-11-27 [1] CRAN (R 4.4.1)
#  lattice                0.22-6    2024-03-20 [1] CRAN (R 4.4.1)
#  lazyeval               0.2.2     2019-03-15 [1] CRAN (R 4.4.0)
#  lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.4.0)
#  limma                  3.60.6    2024-10-02 [1] Bioconductor 3.19 (R 4.4.1)
#  locfit                 1.5-9.10  2024-06-24 [1] CRAN (R 4.4.0)
#  lubridate            * 1.9.4     2024-12-08 [1] CRAN (R 4.4.1)
#  magick                 2.8.5     2024-09-20 [1] CRAN (R 4.4.1)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.4.0)
#  maps                   3.4.2.1   2024-11-10 [1] CRAN (R 4.4.1)
#  Matrix                 1.7-1     2024-10-18 [1] CRAN (R 4.4.1)
#  MatrixGenerics       * 1.16.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  matrixStats          * 1.5.0     2025-01-07 [1] CRAN (R 4.4.1)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.4.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.4.0)
#  munsell                0.5.1     2024-04-01 [1] CRAN (R 4.4.0)
#  paletteer              1.6.0     2024-01-21 [1] CRAN (R 4.4.0)
#  pillar                 1.10.1    2025-01-07 [1] CRAN (R 4.4.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.4.0)
#  plotly                 4.10.4    2024-01-13 [1] CRAN (R 4.4.0)
#  png                    0.1-8     2022-11-29 [1] CRAN (R 4.4.0)
#  promises               1.3.2     2024-11-28 [1] CRAN (R 4.4.1)
#  purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.4.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.4.0)
#  ragg                   1.3.3     2024-09-11 [1] CRAN (R 4.4.1)
#  rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.4.0)
#  RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.4.0)
#  Rcpp                   1.0.14    2025-01-12 [1] CRAN (R 4.4.1)
#  RCurl                  1.98-1.16 2024-07-11 [1] CRAN (R 4.4.0)
#  readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.4.0)
#  rematch2               2.1.2     2020-05-01 [1] CRAN (R 4.4.0)
#  restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.4.0)
#  rjson                  0.2.23    2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.4     2024-06-04 [1] CRAN (R 4.4.0)
#  rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.4.0)
#  Rsamtools              2.20.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  RSQLite                2.3.9     2024-12-03 [1] CRAN (R 4.4.1)
#  rstatix                0.7.2     2023-02-01 [1] CRAN (R 4.4.0)
#  rstudioapi             0.17.1    2024-10-22 [1] CRAN (R 4.4.1)
#  rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.4.0)
#  rtracklayer            1.64.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Arrays               1.4.1     2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Vectors            * 0.42.1    2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  sass                   0.4.9     2024-03-15 [1] CRAN (R 4.4.0)
#  ScaledMatrix           1.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  scales                 1.3.0     2023-11-28 [1] CRAN (R 4.4.0)
#  scater               * 1.32.1    2024-07-21 [1] Bioconductor 3.19 (R 4.4.1)
#  scuttle              * 1.14.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.4.0)
#  shiny                  1.10.0    2024-12-14 [1] CRAN (R 4.4.1)
#  shinyWidgets           0.8.7     2024-09-23 [1] CRAN (R 4.4.1)
#  SingleCellExperiment * 1.26.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  spam                   2.11-0    2024-10-03 [1] CRAN (R 4.4.1)
#  SparseArray            1.4.8     2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sparseMatrixStats      1.16.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SpatialExperiment    * 1.14.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  spatialLIBD          * 1.16.2    2024-05-28 [1] Bioconductor 3.19 (R 4.4.1)
#  statmod                1.5.0     2023-01-06 [1] CRAN (R 4.4.0)
#  stringi                1.8.4     2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.34.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  systemfonts            1.1.0     2024-05-15 [1] CRAN (R 4.4.0)
#  textshaping            0.4.1     2024-12-06 [1] CRAN (R 4.4.1)
#  tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.4.0)
#  tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.4.0)
#  timechange             0.3.0     2024-01-18 [1] CRAN (R 4.4.0)
#  tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.0.0     2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  unigd                  0.1.2     2024-06-05 [1] CRAN (R 4.4.0)
#  vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.4.0)
#  vipor                  0.4.7     2023-12-18 [1] CRAN (R 4.4.0)
#  viridis                0.6.5     2024-01-29 [1] CRAN (R 4.4.0)
#  viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.4.0)
#  withr                  3.0.2     2024-10-28 [1] CRAN (R 4.4.1)
#  XML                    3.99-0.18 2025-01-01 [1] CRAN (R 4.4.1)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.4.0)
#  XVector                0.44.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  yaml                   2.3.10    2024-07-26 [1] CRAN (R 4.4.0)
#  zlibbioc               1.50.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

# ─────────────────────────────────────────────────────────────────────────────────────────────
