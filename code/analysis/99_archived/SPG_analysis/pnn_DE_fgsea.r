# Load library ----
# library(data.table)
library(fgsea)
library(ggplot2)
library(here)
library(tidyverse)
library(sessioninfo)
library(msigdbr)


# Intro to MSigDB ----
# See all collection in MSigbr
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
msigdbr_collections()

## Define MSigDB Gene set ----
msigdbr_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5") # GO genesets
# NOTE: the following allows a more comprehensive gene sets.
# msigdbr_gene_sets <- all_gene_sets |> filter(
#   # gs_cat %in% c("H", "C2", "C3", "C5", "C8")
# )

## Format Gene set for FGSEA ----
# NOTE: we'll use ensembl id for the enrichment analysis
msigdbr_list <- split(
  x = msigdbr_gene_sets$ensembl_gene,
  f = msigdbr_gene_sets$gs_name
)


# Intro to BrainGMT ----
# NOTE: adapted from https://github.com/hagenaue/Brain_GMT/blob/main/BrainGMT_exampleUsage.R
# Load Brain GMT pathway
BrainGMT <- gmtPathways(
  here(
    "code/analysis/SPG_analysis",
    "BrainGMTv2_HumanOrthologs.gmt.txt"
  )
)

# Batch analysis ----
# NOTE: looping over all the DEG csv file
pnn_deg_files <- list.files(here("processed-data/spg_pb_de"))

pnn_deg_files |>
  walk(.f = function(.pb_file) {
    # browser()
    # .pb_file <- "test_SPD_pseudo_pnn_pos.csv"

    ## Load Data ----

    dx_res <- read_csv(
      here(
        "processed-data/spg_pb_de",
        .pb_file
      )
    )

    ## fgsea analysis ----
    ### GO via MSigDB ----
    # Order the genes based on t_statistics in SCZ
    ordered_gene_vector_msigdb <- dx_res |>
      arrange(desc(t_stat_scz)) |>
      select(ensembl, t_stat_scz) |>
      deframe()

    set.seed(20241122)
    msigdb_res <- fgsea(
      pathways = msigdbr_list,
      stats = ordered_gene_vector_msigdb,
      nPermSimple = 10000
    )

    msigdb_res |>
      arrange(padj) |>
      filter(padj < 0.05) |>
      mutate(regulation = case_when(
        NES > 0 ~ "Up in SCZ",
        NES <= 0 ~ "Down in SCZ"
      )) |>
      write_csv(
        here("processed-data/PB_dx_pnn/fgsea", .pb_file)
      )

    ### BrainGMT ----
    set.seed(20241125)
    ordered_gene_vector_braingmt <- dx_res |>
      arrange(desc(t_stat_scz)) |>
      select(gene, t_stat_scz) |>
      deframe()

    braingmt_res <- fgsea(BrainGMT, ordered_gene_vector_braingmt,
      minSize = 10, maxSize = 1000,
      nPermSimple = 10000
    )

    braingmt_res$leadingEdge <- vapply(braingmt_res$leadingEdge,
      paste,
      collapse = ",", character(1L)
    )

    braingmt_res |>
      arrange(padj) |>
      filter(padj < 0.05) |>
      mutate(regulation = case_when(
        NES > 0 ~ "Up in SCZ",
        NES <= 0 ~ "Down in SCZ"
      )) |>
      write_csv(
        here(
          "processed-data/PB_dx_pnn/fgsea",
          str_replace(.pb_file, ".csv", "_brainGMT.csv")
        )
      )
  })


# Session Info ----
sessioninfo::session_info()
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.1 (2024-06-14)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2024-11-22
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────
#  ! package              * version date (UTC) lib source
#    abind                  1.4-8   2024-09-12 [1] CRAN (R 4.4.1)
#    babelgene              22.9    2022-09-29 [1] CRAN (R 4.4.0)
#    beachmat               2.20.0  2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#    beeswarm               0.4.0   2021-06-01 [1] CRAN (R 4.4.0)
#    Biobase              * 2.64.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    BiocGenerics         * 0.50.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    BiocManager            1.30.25 2024-08-28 [1] CRAN (R 4.4.1)
#    BiocNeighbors          1.22.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    BiocParallel           1.38.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    BiocSingular           1.20.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    bit                    4.5.0   2024-09-20 [1] CRAN (R 4.4.1)
#    bit64                  4.5.2   2024-09-22 [1] CRAN (R 4.4.1)
#    cli                    3.6.3   2024-06-21 [1] CRAN (R 4.4.0)
#    codetools              0.2-20  2024-03-31 [1] CRAN (R 4.4.1)
#    colorspace             2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
#    cowplot                1.1.3   2024-01-22 [1] CRAN (R 4.4.0)
#    crayon                 1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
#    data.table           * 1.16.2  2024-10-10 [1] CRAN (R 4.4.1)
#    DelayedArray           0.30.1  2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#    DelayedMatrixStats     1.26.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    dplyr                * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
#    escheR               * 1.4.0   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#    fansi                  1.0.6   2023-12-08 [1] CRAN (R 4.4.0)
#    farver                 2.1.2   2024-05-13 [1] CRAN (R 4.4.0)
#    fastmatch              1.1-4   2023-08-18 [1] CRAN (R 4.4.0)
#    fgsea                * 1.30.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    FNN                    1.1.4.1 2024-09-22 [1] CRAN (R 4.4.1)
#    forcats              * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
#    generics               0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
#    GenomeInfoDb         * 1.40.1  2024-06-16 [1] Bioconductor 3.19 (R 4.4.0)
#    GenomeInfoDbData       1.2.12  2024-08-05 [1] Bioconductor
#  V GenomicRanges        * 1.56.1  2024-10-09 [1] Bioconductor 3.19 (R 4.4.1) (on disk 1.56.2)
#    ggbeeswarm           * 0.7.2   2023-04-29 [1] CRAN (R 4.4.0)
#    ggplot2              * 3.5.1   2024-04-23 [1] CRAN (R 4.4.0)
#    ggrepel                0.9.6   2024-09-07 [1] CRAN (R 4.4.1)
#    glue                   1.8.0   2024-09-30 [1] CRAN (R 4.4.1)
#    gridExtra              2.3     2017-09-09 [1] CRAN (R 4.4.0)
#  V gtable                 0.3.5   2024-10-25 [1] CRAN (R 4.4.1) (on disk 0.3.6)
#    here                 * 1.0.1   2020-12-13 [1] CRAN (R 4.4.0)
#    hms                    1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
#    httpgd                 2.0.2   2024-06-05 [1] CRAN (R 4.4.0)
#    httr                   1.4.7   2023-08-15 [1] CRAN (R 4.4.0)
#    IRanges              * 2.38.1  2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#    irlba                  2.3.5.1 2022-10-03 [1] CRAN (R 4.4.0)
#    jsonlite               1.8.9   2024-09-20 [1] CRAN (R 4.4.1)
#    labeling               0.4.3   2023-08-29 [1] CRAN (R 4.4.0)
#    lattice                0.22-6  2024-03-20 [1] CRAN (R 4.4.1)
#    lifecycle              1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
#    lubridate            * 1.9.3   2023-09-27 [1] CRAN (R 4.4.0)
#  V magick                 2.8.4   2024-09-20 [1] CRAN (R 4.4.1) (on disk 2.8.5)
#    magrittr               2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
#  V Matrix                 1.7-0   2024-10-18 [1] CRAN (R 4.4.1) (on disk 1.7.1)
#    MatrixGenerics       * 1.16.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    matrixStats          * 1.4.1   2024-09-08 [1] CRAN (R 4.4.1)
#    msigdbr              * 7.5.1   2022-03-30 [1] CRAN (R 4.4.0)
#    munsell                0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
#    pillar                 1.9.0   2023-03-22 [1] CRAN (R 4.4.0)
#    pkgconfig              2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
#    purrr                * 1.0.2   2023-08-10 [1] CRAN (R 4.4.0)
#    R6                     2.5.1   2021-08-19 [1] CRAN (R 4.4.0)
#  V Rcpp                   1.0.13  2024-11-02 [1] CRAN (R 4.4.1) (on disk 1.0.13.1)
#    RcppAnnoy              0.0.22  2024-01-23 [1] CRAN (R 4.4.0)
#    readr                * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
#    rjson                  0.2.23  2024-09-16 [1] CRAN (R 4.4.1)
#    rlang                  1.1.4   2024-06-04 [1] CRAN (R 4.4.0)
#    rprojroot              2.0.4   2023-11-05 [1] CRAN (R 4.4.0)
#    rsvd                   1.0.5   2021-04-16 [1] CRAN (R 4.4.0)
#    Rtsne                  0.17    2023-12-07 [1] CRAN (R 4.4.0)
#    S4Arrays               1.4.1   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#    S4Vectors            * 0.42.1  2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
#    ScaledMatrix           1.12.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    scales                 1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
#    scater               * 1.32.1  2024-07-21 [1] Bioconductor 3.19 (R 4.4.1)
#    scuttle              * 1.14.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    sessioninfo          * 1.2.2   2021-12-06 [1] CRAN (R 4.4.0)
#    SingleCellExperiment * 1.26.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    SparseArray            1.4.8   2024-05-30 [1] Bioconductor 3.19 (R 4.4.0)
#    sparseMatrixStats      1.16.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    SpatialExperiment    * 1.14.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    stringi                1.8.4   2024-05-06 [1] CRAN (R 4.4.0)
#    stringr              * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
#    SummarizedExperiment * 1.34.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    systemfonts            1.1.0   2024-05-15 [1] CRAN (R 4.4.0)
#    tibble               * 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
#    tidyr                * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
#    tidyselect             1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
#    tidyverse            * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
#    timechange             0.3.0   2024-01-18 [1] CRAN (R 4.4.0)
#    tzdb                   0.4.0   2023-05-12 [1] CRAN (R 4.4.0)
#    UCSC.utils             1.0.0   2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
#    unigd                  0.1.2   2024-06-05 [1] CRAN (R 4.4.0)
#    utf8                   1.2.4   2023-10-22 [1] CRAN (R 4.4.0)
#    uwot                   0.2.2   2024-04-21 [1] CRAN (R 4.4.0)
#    vctrs                  0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
#    vipor                  0.4.7   2023-12-18 [1] CRAN (R 4.4.0)
#    viridis                0.6.5   2024-01-29 [1] CRAN (R 4.4.0)
#    viridisLite            0.4.2   2023-05-02 [1] CRAN (R 4.4.0)
#    vroom                  1.6.5   2023-12-05 [1] CRAN (R 4.4.0)
#  V withr                  3.0.1   2024-10-28 [1] CRAN (R 4.4.1) (on disk 3.0.2)
#    XVector                0.44.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
#    zlibbioc               1.50.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

#  V ── Loaded and on-disk version mismatch.
