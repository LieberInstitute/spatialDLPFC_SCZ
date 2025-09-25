# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(CellChat)
  library(here)
  library(readxl)
  library(sessioninfo)
})

# Load selected Pathways ----
# Read in Excel file ----
selected_pathways <- read_excel(
  here(
    "code/analysis/Layer_layer_communication",
    "Path_Dx-DEGs.xlsx"
  )
) |> unlist()


db_df <- CellChatDB.human$interaction

ret_df <- db_df |>
  filter(
    pathway_name %in% selected_pathways
  )

ret_df$pathway_name |>
  unique() |>
  length()

write_csv(
  ret_df,
  here(
    "code/analysis/Layer_layer_communication",
    "CellChat_LR_for_Path_Dx-DEGs.csv"
  )
)

# Session info ----
session_info()
# ─ Session info ─────────────────
#  setting  value
#  version  R version 4.4.1 (2024-06-14)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-01-07
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ─────────────────────
#  package              * version  date (UTC) lib source
#  abind                  1.4-8    2024-09-12 [1] CRAN (R 4.4.1)
#  backports              1.5.0    2024-05-23 [1] CRAN (R 4.4.0)
#  Biobase              * 2.64.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocGenerics         * 0.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocManager            1.30.25  2024-08-28 [1] CRAN (R 4.4.1)
#  BiocNeighbors          1.22.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocParallel           1.38.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  bit                    4.5.0    2024-09-20 [1] CRAN (R 4.4.1)
#  bit64                  4.5.2    2024-09-22 [1] CRAN (R 4.4.1)
#  broom                  1.0.7    2024-09-26 [1] CRAN (R 4.4.1)
#  bslib                  0.8.0    2024-07-29 [1] CRAN (R 4.4.0)
#  cachem                 1.1.0    2024-05-16 [1] CRAN (R 4.4.0)
#  car                    3.1-3    2024-09-27 [1] CRAN (R 4.4.1)
#  carData                3.0-5    2022-01-06 [1] CRAN (R 4.4.0)
#  CellChat             * 2.1.2    2024-10-01 [1] Github (jinworks/CellChat@f6aa85d)
#  cellranger             1.1.0    2016-07-27 [1] CRAN (R 4.4.0)
#  circlize               0.4.16   2024-02-20 [1] CRAN (R 4.4.0)
#  cli                    3.6.3    2024-06-21 [1] CRAN (R 4.4.0)
#  clue                   0.3-66   2024-11-13 [1] CRAN (R 4.4.1)
#  cluster                2.1.6    2023-12-01 [1] CRAN (R 4.4.1)
#  coda                   0.19-4.1 2024-01-31 [1] CRAN (R 4.4.0)
#  codetools              0.2-20   2024-03-31 [1] CRAN (R 4.4.1)
#  colorspace             2.1-1    2024-07-26 [1] CRAN (R 4.4.0)
#  ComplexHeatmap         2.20.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  cowplot                1.1.3    2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3    2024-06-20 [1] CRAN (R 4.4.0)
#  DelayedArray           0.30.1   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  digest                 0.6.37   2024-08-19 [1] CRAN (R 4.4.1)
#  doParallel             1.0.17   2022-02-07 [1] CRAN (R 4.4.0)
#  dplyr                * 1.1.4    2023-11-17 [1] CRAN (R 4.4.0)
#  escheR               * 1.4.0    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  fansi                  1.0.6    2023-12-08 [1] CRAN (R 4.4.0)
#  farver                 2.1.2    2024-05-13 [1] CRAN (R 4.4.0)
#  fastmap                1.2.0    2024-05-15 [1] CRAN (R 4.4.0)
#  FNN                    1.1.4.1  2024-09-22 [1] CRAN (R 4.4.1)
#  forcats              * 1.0.0    2023-01-29 [1] CRAN (R 4.4.0)
#  foreach                1.5.2    2022-02-02 [1] CRAN (R 4.4.0)
#  Formula                1.2-5    2023-02-24 [1] CRAN (R 4.4.0)
#  future                 1.34.0   2024-07-29 [1] CRAN (R 4.4.0)
#  future.apply           1.11.3   2024-10-27 [1] CRAN (R 4.4.1)
#  generics               0.1.3    2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.40.1   2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomeInfoDbData       1.2.12   2024-08-05 [1] Bioconductor
#  GenomicRanges        * 1.56.2   2024-10-09 [1] Bioconductor 3.19 (R 4.4.1)
#  GetoptLong             1.0.5    2020-12-15 [1] CRAN (R 4.4.0)
#  ggalluvial             0.12.5   2023-02-22 [1] CRAN (R 4.4.0)
#  ggnetwork              0.5.13   2024-02-14 [1] CRAN (R 4.4.0)
#  ggplot2              * 3.5.1    2024-04-23 [1] CRAN (R 4.4.0)
#  ggpubr               * 0.6.0    2023-02-10 [1] CRAN (R 4.4.0)
#  ggrepel                0.9.6    2024-09-07 [1] CRAN (R 4.4.1)
#  ggsignif               0.6.4    2022-10-13 [1] CRAN (R 4.4.0)
#  GlobalOptions          0.1.2    2020-06-10 [1] CRAN (R 4.4.0)
#  globals                0.16.3   2024-03-08 [1] CRAN (R 4.4.0)
#  glue                   1.8.0    2024-09-30 [1] CRAN (R 4.4.1)
#  gridBase               0.4-7    2014-02-24 [1] CRAN (R 4.4.0)
#  gridExtra              2.3      2017-09-09 [1] CRAN (R 4.4.0)
#  gtable                 0.3.6    2024-10-25 [1] CRAN (R 4.4.1)
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3    2023-03-21 [1] CRAN (R 4.4.0)
#  htmltools              0.5.8.1  2024-04-04 [1] CRAN (R 4.4.0)
#  httpuv                 1.6.15   2024-03-26 [1] CRAN (R 4.4.0)
#  httr                   1.4.7    2023-08-15 [1] CRAN (R 4.4.0)
#  igraph               * 2.1.1    2024-10-19 [1] CRAN (R 4.4.1)
#  IRanges              * 2.38.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  irlba                  2.3.5.1  2022-10-03 [1] CRAN (R 4.4.0)
#  iterators              1.0.14   2022-02-05 [1] CRAN (R 4.4.0)
#  jquerylib              0.1.4    2021-04-26 [1] CRAN (R 4.4.0)
#  jsonlite               1.8.9    2024-09-20 [1] CRAN (R 4.4.1)
#  labeling               0.4.3    2023-08-29 [1] CRAN (R 4.4.0)
#  later                  1.3.2    2023-12-06 [1] CRAN (R 4.4.0)
#  lattice                0.22-6   2024-03-20 [1] CRAN (R 4.4.1)
#  lifecycle              1.0.4    2023-11-07 [1] CRAN (R 4.4.0)
#  listenv                0.9.1    2024-01-29 [1] CRAN (R 4.4.0)
#  lubridate            * 1.9.3    2023-09-27 [1] CRAN (R 4.4.0)
#  magick                 2.8.5    2024-09-20 [1] CRAN (R 4.4.1)
#  magrittr               2.0.3    2022-03-30 [1] CRAN (R 4.4.0)
#  Matrix                 1.7-1    2024-10-18 [1] CRAN (R 4.4.1)
#  MatrixGenerics       * 1.16.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  matrixStats          * 1.4.1    2024-09-08 [1] CRAN (R 4.4.1)
#  mime                   0.12     2021-09-28 [1] CRAN (R 4.4.0)
#  munsell                0.5.1    2024-04-01 [1] CRAN (R 4.4.0)
#  network                1.18.2   2023-12-05 [1] CRAN (R 4.4.0)
#  NMF                    0.28     2024-08-22 [1] CRAN (R 4.4.1)
#  parallelly             1.39.0   2024-11-07 [1] CRAN (R 4.4.1)
#  patchwork              1.3.0    2024-09-16 [1] CRAN (R 4.4.1)
#  pbapply                1.7-2    2023-06-27 [1] CRAN (R 4.4.0)
#  pillar                 1.9.0    2023-03-22 [1] CRAN (R 4.4.0)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.4.0)
#  plyr                   1.8.9    2023-10-02 [1] CRAN (R 4.4.0)
#  png                    0.1-8    2022-11-29 [1] CRAN (R 4.4.0)
#  Polychrome             1.5.1    2022-05-03 [1] CRAN (R 4.4.0)
#  promises               1.3.0    2024-04-05 [1] CRAN (R 4.4.0)
#  purrr                * 1.0.2    2023-08-10 [1] CRAN (R 4.4.0)
#  R6                     2.5.1    2021-08-19 [1] CRAN (R 4.4.0)
#  RColorBrewer           1.1-3    2022-04-03 [1] CRAN (R 4.4.0)
#  Rcpp                   1.0.13-1 2024-11-02 [1] CRAN (R 4.4.1)
#  readr                * 2.1.5    2024-01-10 [1] CRAN (R 4.4.0)
#  readxl               * 1.4.3    2023-07-06 [1] CRAN (R 4.4.0)
#  registry               0.5-1    2019-03-05 [1] CRAN (R 4.4.0)
#  reshape2               1.4.4    2020-04-09 [1] CRAN (R 4.4.0)
#  reticulate             1.40.0   2024-11-15 [1] CRAN (R 4.4.1)
#  rjson                  0.2.23   2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.4    2024-06-04 [1] CRAN (R 4.4.0)
#  rngtools               1.5.2    2021-09-20 [1] CRAN (R 4.4.0)
#  rprojroot              2.0.4    2023-11-05 [1] CRAN (R 4.4.0)
#  RSpectra               0.16-2   2024-07-18 [1] CRAN (R 4.4.0)
#  rstatix                0.7.2    2023-02-01 [1] CRAN (R 4.4.0)
#  S4Arrays               1.4.1    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Vectors            * 0.42.1   2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  sass                   0.4.9    2024-03-15 [1] CRAN (R 4.4.0)
#  scales                 1.3.0    2023-11-28 [1] CRAN (R 4.4.0)
#  scatterplot3d          0.3-44   2023-05-05 [1] CRAN (R 4.4.0)
#  sessioninfo          * 1.2.2    2021-12-06 [1] CRAN (R 4.4.0)
#  shape                  1.4.6.1  2024-02-23 [1] CRAN (R 4.4.0)
#  shiny                  1.9.1    2024-08-01 [1] CRAN (R 4.4.0)
#  SingleCellExperiment * 1.26.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sna                    2.8      2024-09-08 [1] CRAN (R 4.4.1)
#  SparseArray            1.4.8    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SpatialExperiment    * 1.14.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  statnet.common         4.10.0   2024-10-06 [1] CRAN (R 4.4.1)
#  stringi                1.8.4    2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1    2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.34.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  svglite                2.1.3    2023-12-08 [1] CRAN (R 4.4.0)
#  systemfonts            1.1.0    2024-05-15 [1] CRAN (R 4.4.0)
#  tibble               * 3.2.1    2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                * 1.3.1    2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1    2024-03-11 [1] CRAN (R 4.4.0)
#  tidyverse            * 2.0.0    2023-02-22 [1] CRAN (R 4.4.0)
#  timechange             0.3.0    2024-01-18 [1] CRAN (R 4.4.0)
#  tzdb                   0.4.0    2023-05-12 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.0.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  utf8                   1.2.4    2023-10-22 [1] CRAN (R 4.4.0)
#  vctrs                  0.6.5    2023-12-01 [1] CRAN (R 4.4.0)
#  vroom                  1.6.5    2023-12-05 [1] CRAN (R 4.4.0)
#  withr                  3.0.2    2024-10-28 [1] CRAN (R 4.4.1)
#  xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.4.0)
#  XVector                0.44.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  zlibbioc               1.50.0   2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library