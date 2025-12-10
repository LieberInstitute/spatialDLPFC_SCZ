# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
  library(ggrepel)
})

# Load interaction fit ----
fit <- readRDS(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "limma_obj_int_PRECAST_07.rds"
  )
)


# Create contrast matrix of SpDs ----
cont.mat <- rbind(
  rep(-1, 7),
  rep(1, 7),
  matrix(0, nrow = 23, ncol = 7),
  cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
)
colnames(cont.mat) <- sprintf("spd%02d", 1:7)

## Output of contrasts ----
#       spd01 spd02 spd03 spd04 spd05 spd06 spd07
#  [1,]    -1    -1    -1    -1    -1    -1    -1
#  [2,]     1     1     1     1     1     1     1
#  [3,]     0     0     0     0     0     0     0
#  [4,]     0     0     0     0     0     0     0
#  [5,]     0     0     0     0     0     0     0
#  [6,]     0     0     0     0     0     0     0
#  [7,]     0     0     0     0     0     0     0
#  [8,]     0     0     0     0     0     0     0
#  [9,]     0     0     0     0     0     0     0
# [10,]     0     0     0     0     0     0     0
# [11,]     0     0     0     0     0     0     0
# [12,]     0     0     0     0     0     0     0
# [13,]     0     0     0     0     0     0     0
# [14,]     0     0     0     0     0     0     0
# [15,]     0     0     0     0     0     0     0
# [16,]     0     0     0     0     0     0     0
# [17,]     0     0     0     0     0     0     0
# [18,]     0     0     0     0     0     0     0
# [19,]     0     0     0     0     0     0     0
# [20,]     0     0     0     0     0     0     0
# [21,]     0     0     0     0     0     0     0
# [22,]     0     0     0     0     0     0     0
# [23,]     0     0     0     0     0     0     0
# [24,]     0     0     0     0     0     0     0
# [25,]     0     0     0     0     0     0     0
# [26,]     0     1     0     0     0     0     0
# [27,]     0     0     1     0     0     0     0
# [28,]     0     0     0     1     0     0     0
# [29,]     0     0     0     0     1     0     0
# [30,]     0     0     0     0     0     1     0
# [31,]     0     0     0     0     0     0     1

# Error prevention
# Design matrix is conformable with the contrast matrix
stopifnot(nrow(cont.mat) == ncol(fit$design))

# Layer specific logFC and hypo test ----

colnames(cont.mat) |>
  set_names() |>
  iwalk(.f = function(spd_it, idx) {
    # browser()
    spd_cont_df <- topTable(
      contrast_fit,
      coef = spd_it, num = Inf
    ) |>
      rownames_to_column("gene_id") |>
      mutate(
        sig_10 = adj.P.Val <= 0.1
      )

    write_csv(
      spd_cont_df,
      here(
        "processed-data/rds/11_dx_deg_interaction",
        paste0(
          # "layer_specific_logFC_",
          "layer_restricted_logFC_",
          gsub("[^[:alnum:]_]", "_", idx),
          ".csv"
        )
      )
    )
  })




# Session Info ----
# sessioninfo::session_info()
# ─ Session info ─────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.2 (2024-10-31)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-04-22
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ─────────────────────────────────────────────────
#  package              * version date (UTC) lib source
#  abind                  1.4-8   2024-09-12 [1] CRAN (R 4.4.1)
#  AnnotationDbi          1.68.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  ape                    5.8-1   2024-12-16 [1] CRAN (R 4.4.1)
#  aplot                  0.2.5   2025-02-27 [1] CRAN (R 4.4.1)
#  beachmat               2.22.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  beeswarm               0.4.0   2021-06-01 [1] CRAN (R 4.4.0)
#  Biobase              * 2.66.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocGenerics         * 0.52.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocNeighbors          2.0.1   2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
#  BiocParallel           1.40.0  2024-11-23 [1] Bioconductor 3.20 (R 4.4.2)
#  BiocSingular           1.22.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  Biostrings             2.74.1  2024-12-16 [1] Bioconductor 3.20 (R 4.4.2)
#  bit                    4.5.0.1 2024-12-03 [1] CRAN (R 4.4.1)
#  bit64                  4.6.0-1 2025-01-16 [1] CRAN (R 4.4.1)
#  blob                   1.2.4   2023-03-17 [1] CRAN (R 4.4.0)
#  cachem                 1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
#  circlize               0.4.16  2024-02-20 [1] CRAN (R 4.4.0)
#  cli                    3.6.3   2024-06-21 [1] CRAN (R 4.4.0)
#  clue                   0.3-66  2024-11-13 [1] CRAN (R 4.4.1)
#  cluster                2.1.6   2023-12-01 [1] CRAN (R 4.4.2)
#  clusterProfiler        4.14.6  2025-02-27 [1] Bioconductor 3.20 (R 4.4.2)
#  codetools              0.2-20  2024-03-31 [1] CRAN (R 4.4.2)
#  colorspace             2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
#  ComplexHeatmap       * 2.22.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  cowplot                1.1.3   2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
#  data.table             1.16.4  2024-12-06 [1] CRAN (R 4.4.1)
#  DBI                    1.2.3   2024-06-02 [1] CRAN (R 4.4.0)
#  DelayedArray           0.32.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  digest                 0.6.37  2024-08-19 [1] CRAN (R 4.4.1)
#  doParallel             1.0.17  2022-02-07 [1] CRAN (R 4.4.0)
#  DOSE                   4.0.1   2025-03-27 [1] Bioconductor 3.20 (R 4.4.3)
#  dplyr                * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
#  enrichplot             1.26.6  2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
#  farver                 2.1.2   2024-05-13 [1] CRAN (R 4.4.0)
#  fastmap                1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
#  fastmatch              1.1-6   2024-12-23 [1] CRAN (R 4.4.1)
#  fgsea                  1.32.2  2024-12-19 [1] Bioconductor 3.20 (R 4.4.1)
#  forcats              * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
#  foreach                1.5.2   2022-02-02 [1] CRAN (R 4.4.0)
#  fs                     1.6.5   2024-10-30 [1] CRAN (R 4.4.1)
#  generics               0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.42.1  2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
#  GenomeInfoDbData       1.2.13  2025-01-19 [1] Bioconductor
#  GenomicRanges        * 1.58.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  GetoptLong             1.0.5   2020-12-15 [1] CRAN (R 4.4.0)
#  ggbeeswarm             0.7.2   2023-04-29 [1] CRAN (R 4.4.0)
#  ggfun                  0.1.8   2024-12-03 [1] CRAN (R 4.4.1)
#  ggplot2              * 3.5.1   2024-04-23 [1] CRAN (R 4.4.0)
#  ggplotify              0.1.2   2023-08-09 [1] CRAN (R 4.4.0)
#  ggrepel              * 0.9.6   2024-09-07 [1] CRAN (R 4.4.1)
#  ggtangle               0.0.6   2024-12-18 [1] CRAN (R 4.4.1)
#  ggtree                 3.14.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  GlobalOptions          0.1.2   2020-06-10 [1] CRAN (R 4.4.0)
#  glue                   1.8.0   2024-09-30 [1] CRAN (R 4.4.1)
#  GO.db                  3.20.0  2025-04-16 [1] Bioconductor
#  GOSemSim               2.32.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  gridExtra              2.3     2017-09-09 [1] CRAN (R 4.4.0)
#  gridGraphics           0.5-1   2020-12-13 [1] CRAN (R 4.4.1)
#  gson                   0.1.0   2023-03-07 [1] CRAN (R 4.4.0)
#  gtable                 0.3.6   2024-10-25 [1] CRAN (R 4.4.1)
#  here                 * 1.0.1   2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
#  httr                   1.4.7   2023-08-15 [1] CRAN (R 4.4.0)
#  igraph                 2.1.3   2025-01-07 [1] CRAN (R 4.4.1)
#  IRanges              * 2.40.1  2024-12-05 [1] Bioconductor 3.20 (R 4.4.2)
#  irlba                  2.3.5.1 2022-10-03 [1] CRAN (R 4.4.0)
#  iterators              1.0.14  2022-02-05 [1] CRAN (R 4.4.0)
#  jsonlite               1.8.9   2024-09-20 [1] CRAN (R 4.4.1)
#  KEGGREST               1.46.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  lattice                0.22-6  2024-03-20 [1] CRAN (R 4.4.2)
#  lazyeval               0.2.2   2019-03-15 [1] CRAN (R 4.4.0)
#  lifecycle              1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
#  limma                * 3.62.2  2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
#  lubridate            * 1.9.4   2024-12-08 [1] CRAN (R 4.4.1)
#  magick                 2.8.5   2024-09-20 [1] CRAN (R 4.4.1)
#  magrittr               2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
#  Matrix                 1.7-1   2024-10-18 [1] CRAN (R 4.4.2)
#  MatrixGenerics       * 1.18.1  2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
#  matrixStats          * 1.5.0   2025-01-07 [1] CRAN (R 4.4.1)
#  memoise                2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
#  munsell                0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
#  nlme                   3.1-166 2024-08-14 [1] CRAN (R 4.4.2)
#  org.Hs.eg.db           3.20.0  2025-04-16 [1] Bioconductor
#  patchwork              1.3.0   2024-09-16 [1] CRAN (R 4.4.1)
#  pillar                 1.10.1  2025-01-07 [1] CRAN (R 4.4.1)
#  pkgconfig              2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
#  plyr                   1.8.9   2023-10-02 [1] CRAN (R 4.4.0)
#  png                    0.1-8   2022-11-29 [1] CRAN (R 4.4.0)
#  purrr                * 1.0.2   2023-08-10 [1] CRAN (R 4.4.0)
#  qvalue                 2.38.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  R.methodsS3            1.8.2   2022-06-13 [1] CRAN (R 4.4.0)
#  R.oo                   1.27.0  2024-11-01 [1] CRAN (R 4.4.1)
#  R.utils                2.12.3  2023-11-18 [1] CRAN (R 4.4.0)
#  R6                     2.5.1   2021-08-19 [1] CRAN (R 4.4.0)
#  RColorBrewer           1.1-3   2022-04-03 [1] CRAN (R 4.4.0)
#  Rcpp                   1.0.14  2025-01-12 [1] CRAN (R 4.4.1)
#  readr                * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
#  reshape2               1.4.4   2020-04-09 [1] CRAN (R 4.4.0)
#  rjson                  0.2.23  2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.5   2025-01-17 [1] CRAN (R 4.4.1)
#  rprojroot              2.0.4   2023-11-05 [1] CRAN (R 4.4.0)
#  RSQLite                2.3.9   2024-12-03 [1] CRAN (R 4.4.1)
#  rsvd                   1.0.5   2021-04-16 [1] CRAN (R 4.4.0)
#  S4Arrays               1.6.0   2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  S4Vectors            * 0.44.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  ScaledMatrix           1.14.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  scales                 1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
#  scater               * 1.34.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
#  scuttle              * 1.16.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  sessioninfo          * 1.2.2   2021-12-06 [1] CRAN (R 4.4.0)
#  shape                  1.4.6.1 2024-02-23 [1] CRAN (R 4.4.0)
#  SingleCellExperiment * 1.28.1  2024-11-23 [1] Bioconductor 3.20 (R 4.4.2)
#  SparseArray            1.6.0   2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  SpatialExperiment    * 1.16.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  statmod                1.5.0   2023-01-06 [1] CRAN (R 4.4.0)
#  stringi                1.8.4   2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.36.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  tibble               * 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
#  tidytree               0.4.6   2023-12-12 [1] CRAN (R 4.4.0)
#  tidyverse            * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
#  timechange             0.3.0   2024-01-18 [1] CRAN (R 4.4.0)
#  treeio                 1.30.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  tzdb                   0.4.0   2023-05-12 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.2.0   2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  vctrs                  0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
#  vipor                  0.4.7   2023-12-18 [1] CRAN (R 4.4.0)
#  viridis                0.6.5   2024-01-29 [1] CRAN (R 4.4.0)
#  viridisLite            0.4.2   2023-05-02 [1] CRAN (R 4.4.0)
#  vroom                  1.6.5   2023-12-05 [1] CRAN (R 4.4.0)
#  withr                  3.0.2   2024-10-28 [1] CRAN (R 4.4.1)
#  XVector                0.46.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  yulab.utils            0.2.0   2025-01-29 [1] CRAN (R 4.4.1)
#  zlibbioc               1.52.0  2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

# ───────
