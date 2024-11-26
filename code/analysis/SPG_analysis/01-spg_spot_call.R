# Load packages----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

## Load SpD data ----
finalized_spd <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Attach SpD label to spe ----
col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)


# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)

# Call PNN & PVALB + ----
pvalb_index <- which(rowData(spe)$gene_name == "PVALB")
spe$pvalb_pos <- logcounts(spe)[pvalb_index, ] > 0
spe$logcount_pvalb <- logcounts(spe)[pvalb_index, ]

spe$pnn_pvalb <- spe$pnn_pos & spe$pvalb_pos

## Call SPG Neighbors ----
### Prepare data for BayesSpace ----
# In order to use BayesSpace to identify neighbors,
# we first need to have `row` and `col` in the colData
spe$row <- spe$array_row
spe$col <- spe$array_col

### Helpfer function ----
# NOTE: extract neighbors of specific SPG channel
which_neighbors <- function(spe, neighbors_list, var, return_keys = TRUE) {
  i <- which(colData(spe)[[var]] == TRUE)

  # NOTE: the trailing +1 is very important. Don't change
  res <- sort(unique(unlist(neighbors_list[i]))) + 1
  res_keys <- spe$key[res]


  if (return_keys == TRUE) {
    return(res_keys)
  } else {
    return(res)
  }
}

## Calculate PNN Neighbors ----

find_neighbors_for <- function(spe, var) {
  neighbor_key <- c()
  stopifnot( var %in% names(colData(spe)))
  for (.sample in unique(spe$sample_id)) {
    sub_spe <- spe[, spe$sample_id == .sample]
    neighbors_list <- BayesSpace:::.find_neighbors(
      sub_spe,
      platform = "Visium"
    )
    # browser()
    sub_key <- which_neighbors(sub_spe, neighbors_list, var, return_keys = TRUE)
    neighbor_key <- c(neighbor_key, sub_key) # concatenate keys
  }

  neighbor_vec <- rep(FALSE, ncol(spe))
  neighbor_vec[spe$key %in% neighbor_key] <- TRUE

  return(spe[[var]] | neighbor_vec)
}

spe$pnn_N_neighbors <- find_neighbors_for(
  spe = spe, var = "pnn_pos"
)

# Calculate PNN PVALB Neighbors ----
spe$pnn_pvalb_N_neighbors <- find_neighbors_for(
  spe = spe, var = "pnn_pvalb"
)



# Create channel specific PB ----

# spg_names <- c("pnn_pos", "neuropil_pos", "neun_pos", "vasc_pos")
# spg_names <- c("pnn_N_neighbors")
spg_names <- c("pnn_pvalb", "pnn_pvalb_N_neighbors")

## Create positive PNN only spe ----
for (.spg in spg_names) {
  spg_spe <- spe[, spe[[.spg]] == TRUE]

  ## Pseudobulk based on individual
  sce_pseudo <-
    registration_pseudobulk(
      spg_spe,
      var_registration = "dx",
      var_sample_id = "sample_id",
      covars = c("age", "sex", "lot_num", "slide_id"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "rds", "PB_dx_spg",
        paste0("test_donor_pseudo_", .spg, ".rds")
      )
    )

  ## Per sample & domain ----
  sce_pseudo <-
    registration_pseudobulk(
      spg_spe,
      var_registration = "PRECAST_07",
      var_sample_id = "sample_id",
      covars = c("dx", "age", "sex", "lot_num", "slide_id"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "rds", "PB_dx_spg",
        paste0("test_SPD_pseudo_", .spg, ".rds")
      )
    )
}


# Session Info ----
sessioninfo::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.1 (2024-06-14)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2024-09-19
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ───────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [1] CRAN (R 4.4.0)
#  AnnotationDbi          1.66.0    2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  AnnotationHub          3.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  attempt                0.3.1     2020-05-03 [1] CRAN (R 4.4.0)
#  beachmat               2.20.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.4.0)
#  benchmarkme            1.0.8     2022-06-12 [1] CRAN (R 4.4.0)
#  benchmarkmeData        1.0.4     2020-04-23 [1] CRAN (R 4.4.0)
#  Biobase              * 2.64.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocFileCache          2.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocGenerics         * 0.50.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocIO                 1.14.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocManager            1.30.23   2024-05-04 [1] CRAN (R 4.4.0)
#  BiocNeighbors          1.22.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocParallel           1.38.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocSingular           1.20.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  BiocVersion            3.19.1    2024-04-22 [1] Bioconductor 3.19 (R 4.4.0)
#  Biostrings             2.72.1    2024-06-02 [1] Bioconductor 3.19 (R 4.4.0)
#  bit                    4.0.5     2022-11-15 [1] CRAN (R 4.4.0)
#  bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.4.0)
#  bitops                 1.0-8     2024-07-29 [1] CRAN (R 4.4.0)
#  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.4.0)
#  bslib                  0.8.0     2024-07-29 [1] CRAN (R 4.4.0)
#  cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.4.0)
#  cli                    3.6.3     2024-06-21 [1] CRAN (R 4.4.0)
#  codetools              0.2-20    2024-03-31 [1] CRAN (R 4.4.1)
#  colorspace             2.1-1     2024-07-26 [1] CRAN (R 4.4.0)
#  config                 0.3.2     2023-08-30 [1] CRAN (R 4.4.0)
#  cowplot                1.1.3     2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3     2024-06-20 [1] CRAN (R 4.4.0)
#  curl                   5.2.1     2024-03-01 [1] CRAN (R 4.4.0)
#  data.table             1.15.4    2024-03-30 [1] CRAN (R 4.4.0)
#  DBI                    1.2.3     2024-06-02 [1] CRAN (R 4.4.0)
#  dbplyr                 2.5.0     2024-03-19 [1] CRAN (R 4.4.0)
#  DelayedArray           0.30.1    2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  DelayedMatrixStats     1.26.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  digest                 0.6.36    2024-06-23 [1] CRAN (R 4.4.0)
#  doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.4.0)
#  dotCall64              1.1-1     2023-11-28 [1] CRAN (R 4.4.0)
#  dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.4.0)
#  DT                     0.33      2024-04-04 [1] CRAN (R 4.4.0)
#  edgeR                  4.2.1     2024-07-14 [1] Bioconductor 3.19 (R 4.4.1)
#  escheR               * 1.4.0     2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  ExperimentHub          2.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.4.0)
#  fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.4.0)
#  fields                 16.2      2024-06-27 [1] CRAN (R 4.4.0)
#  filelock               1.0.3     2023-12-11 [1] CRAN (R 4.4.0)
#  forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.4.0)
#  foreach                1.5.2     2022-02-02 [1] CRAN (R 4.4.0)
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.40.1    2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomeInfoDbData       1.2.12    2024-08-05 [1] Bioconductor
#  GenomicAlignments      1.40.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  GenomicRanges        * 1.56.1    2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#  ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.4.0)
#  ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.4.0)
#  ggrepel                0.9.5     2024-01-10 [1] CRAN (R 4.4.0)
#  glue                   1.7.0     2024-01-09 [1] CRAN (R 4.4.0)
#  golem                  0.4.1     2023-06-05 [1] CRAN (R 4.4.0)
#  gridExtra              2.3       2017-09-09 [1] CRAN (R 4.4.0)
#  gtable                 0.3.5     2024-04-22 [1] CRAN (R 4.4.0)
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.4.0)
#  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.4.0)
#  htmltools              0.5.8.1   2024-04-04 [1] CRAN (R 4.4.0)
#  htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.4.0)
#  httpuv                 1.6.15    2024-03-26 [1] CRAN (R 4.4.0)
#  httr                   1.4.7     2023-08-15 [1] CRAN (R 4.4.0)
#  IRanges              * 2.38.1    2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.4.0)
#  iterators              1.0.14    2022-02-05 [1] CRAN (R 4.4.0)
#  jquerylib              0.1.4     2021-04-26 [1] CRAN (R 4.4.0)
#  jsonlite               1.8.8     2023-12-04 [1] CRAN (R 4.4.0)
#  KEGGREST               1.44.1    2024-06-19 [1] Bioconductor 3.19 (R 4.4.0)
#  later                  1.3.2     2023-12-06 [1] CRAN (R 4.4.0)
#  lattice                0.22-6    2024-03-20 [1] CRAN (R 4.4.1)
#  lazyeval               0.2.2     2019-03-15 [1] CRAN (R 4.4.0)
#  lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.4.0)
#  limma                * 3.60.4    2024-07-17 [1] Bioconductor 3.19 (R 4.4.1)
#  locfit                 1.5-9.10  2024-06-24 [1] CRAN (R 4.4.0)
#  lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.4.0)
#  magick                 2.8.4     2024-07-14 [1] CRAN (R 4.4.0)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.4.0)
#  maps                   3.4.2     2023-12-15 [1] CRAN (R 4.4.0)
#  Matrix                 1.7-0     2024-04-26 [1] CRAN (R 4.4.1)
#  MatrixGenerics       * 1.16.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  matrixStats          * 1.3.0     2024-04-11 [1] CRAN (R 4.4.0)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.4.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.4.0)
#  munsell                0.5.1     2024-04-01 [1] CRAN (R 4.4.0)
#  paletteer              1.6.0     2024-01-21 [1] CRAN (R 4.4.0)
#  pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.4.0)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.4.0)
#  plotly                 4.10.4    2024-01-13 [1] CRAN (R 4.4.0)
#  png                    0.1-8     2022-11-29 [1] CRAN (R 4.4.0)
#  promises               1.3.0     2024-04-05 [1] CRAN (R 4.4.0)
#  purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.4.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.4.0)
#  rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.4.0)
#  RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.4.0)
#  Rcpp                   1.0.13    2024-07-17 [1] CRAN (R 4.4.0)
#  RCurl                  1.98-1.16 2024-07-11 [1] CRAN (R 4.4.0)
#  readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.4.0)
#  rematch2               2.1.2     2020-05-01 [1] CRAN (R 4.4.0)
#  restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.4.0)
#  rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.4.0)
#  rlang                  1.1.4     2024-06-04 [1] CRAN (R 4.4.0)
#  rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.4.0)
#  Rsamtools              2.20.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  RSQLite                2.3.7     2024-05-27 [1] CRAN (R 4.4.0)
#  rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.4.0)
#  rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.4.0)
#  rtracklayer            1.64.0    2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Arrays               1.4.1     2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  S4Vectors            * 0.42.1    2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#  sass                   0.4.9     2024-03-15 [1] CRAN (R 4.4.0)
#  ScaledMatrix           1.12.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  scales                 1.3.0     2023-11-28 [1] CRAN (R 4.4.0)
#  scater                 1.32.1    2024-07-21 [1] Bioconductor 3.19 (R 4.4.1)
#  scuttle                1.14.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.4.0)
#  shiny                  1.9.1     2024-08-01 [1] CRAN (R 4.4.0)
#  shinyWidgets           0.8.6     2024-04-24 [1] CRAN (R 4.4.0)
#  SingleCellExperiment * 1.26.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  spam                   2.10-0    2023-10-23 [1] CRAN (R 4.4.0)
#  SparseArray            1.4.8     2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#  sparseMatrixStats      1.16.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  SpatialExperiment    * 1.14.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  spatialLIBD          * 1.16.2    2024-05-28 [1] Bioconductor 3.19 (R 4.4.1)
#  statmod                1.5.0     2023-01-06 [1] CRAN (R 4.4.0)
#  stringi                1.8.4     2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.34.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.4.0)
#  tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.4.0)
#  timechange             0.3.0     2024-01-18 [1] CRAN (R 4.4.0)
#  tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.0.0     2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#  utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.4.0)
#  vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.4.0)
#  vipor                  0.4.7     2023-12-18 [1] CRAN (R 4.4.0)
#  viridis                0.6.5     2024-01-29 [1] CRAN (R 4.4.0)
#  viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.4.0)
#  withr                  3.0.1     2024-07-31 [1] CRAN (R 4.4.0)
#  XML                    3.99-0.17 2024-06-25 [1] CRAN (R 4.4.0)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.4.0)
#  XVector                0.44.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#  yaml                   2.3.10    2024-07-26 [1] CRAN (R 4.4.0)
#  zlibbioc               1.50.0    2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

# ──────────────────────────────────────────────────────────────────────────────────
