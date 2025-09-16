suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(sessioninfo)
})

# Load Data ----
## Load SCZ-DEG PB results ----
gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd/preliminary",
    "test_PRECAST_07.csv"
  )
  # here(
  #   "processed-data/rds/10_dx_deg_adjust_spd",
  #   TODO: correct the file name
  #   "test_PRECAST_07.csv"
  # )
) |> mutate(
  sig_05 = fdr_scz <= 0.05,
  sig_10 = fdr_scz <= 0.10,
  gene_cat = case_when(
    sig_10 == FALSE ~ "grey",
    sig_10 == TRUE & logFC_scz >= 0 ~ "red",
    sig_10 == TRUE & logFC_scz < 0 ~ "blue"
  )
)

## Load interesteded gene list ----
gene_names <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    "Escher_genelist.xlsx"
  ),
  col_names = FALSE
) |> unlist()
names(gene_names) <- "gene"

# Code for exploratory plot making
# sig_gene_df <- gene_df |> filter(gene %in% gene_names)

# Gene curated by Sang Ho
sig_gene_df <- gene_df |> filter(gene %in% c("BDNF", "FKBP5", "SST", "GFAP", "MT2A"))

# Make volcano plot ----
ggplot(
  data = gene_df |> arrange(desc(p_value_scz)),
  aes(x = logFC_scz, y = -log10(p_value_scz))
) +
  # Add line for nominal threshold 0.05
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "darkgrey"
  ) +
  # Plot all genes
  geom_point(
    aes(color = gene_cat),
    size = 0.001
  ) +
  # Add interested gene labels
  geom_text_repel(
    data = sig_gene_df, # Add labels last to appear as the top layer
    aes(
      x = logFC_scz, y = -log10(p_value_scz), label = gene # , color = gene_cat
    ),
    color = "black",
    force = 0.5,
    # nudge_y = 0.1,
    # All labels not overlapping
    max.overlaps = Inf,
    # Have arrows
    min.segment.length = 0,
    size = 1.5, # 6pt label font ≈ size 2 in ggplot2
    arrow = arrow(length = unit(0.5, "points"), type = "closed"),
    segment.size = 0.2,
    segment.color = "black"
  ) +
  # Format the legends
  scale_color_identity(
    name = "Significance",
    labels = c(
      "FDR > 0.10",
      "Up-regulated",
      "Down-regulated"
    )
  ) +
  labs(
    title = "Layer-adjusted DEGs",
    y = "-log10(p-value)",
    x = "log2(FC in SCZ)",
    color = "Significance"
  ) +
  # Format the plot
  theme_classic(base_size = 6) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    #   axis.title = element_text(size = 6),
    #   axis.text = element_text(size = 6),
    #   legend.position = "top",
    #   panel.border = element_rect(
    #     color = "black", fill = NA, size = 1
    #   )
  )
#
#
# Save the plot ----
ggsave(
  filename = here(
    "plots/10_dx_deg_adjust_spd",
    "volcano_PRECAST_07.pdf"
  ),
  width = 1.72,
  height = 1.3,
  units = "in",
  dpi = 300
)

# Session Info ----
session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.4.2 (2024-10-31)
#  os       macOS Sonoma 14.6.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2025-03-19
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc

# ─ Packages ──────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-8     2024-09-12 [1] CRAN (R 4.4.1)
#  AnnotationDbi          1.68.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  AnnotationHub          3.14.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  attempt                0.3.1     2020-05-03 [1] CRAN (R 4.4.0)
#  beachmat               2.22.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.4.0)
#  benchmarkme            1.0.8     2022-06-12 [1] CRAN (R 4.4.0)
#  benchmarkmeData        1.0.4     2020-04-23 [1] CRAN (R 4.4.0)
#  Biobase              * 2.66.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocFileCache          2.14.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocGenerics         * 0.52.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocIO                 1.16.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocManager            1.30.25   2024-08-28 [1] CRAN (R 4.4.1)
#  BiocNeighbors          2.0.1     2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
#  BiocParallel           1.40.0    2024-11-23 [1] Bioconductor 3.20 (R 4.4.2)
#  BiocSingular           1.22.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  BiocVersion            3.20.0    2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
#  Biostrings             2.74.1    2024-12-16 [1] Bioconductor 3.20 (R 4.4.2)
#  bit                    4.5.0.1   2024-12-03 [1] CRAN (R 4.4.1)
#  bit64                  4.6.0-1   2025-01-16 [1] CRAN (R 4.4.1)
#  bitops                 1.0-9     2024-10-03 [1] CRAN (R 4.4.1)
#  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.4.0)
#  bslib                  0.8.0     2024-07-29 [1] CRAN (R 4.4.0)
#  cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.4.0)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.4.0)
#  cli                    3.6.3     2024-06-21 [1] CRAN (R 4.4.0)
#  codetools              0.2-20    2024-03-31 [1] CRAN (R 4.4.2)
#  colorspace             2.1-1     2024-07-26 [1] CRAN (R 4.4.0)
#  config                 0.3.2     2023-08-30 [1] CRAN (R 4.4.0)
#  cowplot                1.1.3     2024-01-22 [1] CRAN (R 4.4.0)
#  crayon                 1.5.3     2024-06-20 [1] CRAN (R 4.4.0)
#  curl                   6.1.0     2025-01-06 [1] CRAN (R 4.4.2)
#  data.table             1.16.4    2024-12-06 [1] CRAN (R 4.4.1)
#  DBI                    1.2.3     2024-06-02 [1] CRAN (R 4.4.0)
#  dbplyr                 2.5.0     2024-03-19 [1] CRAN (R 4.4.0)
#  DelayedArray           0.32.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  digest                 0.6.37    2024-08-19 [1] CRAN (R 4.4.1)
#  doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.4.0)
#  dotCall64              1.2       2024-10-04 [1] CRAN (R 4.4.1)
#  dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.4.0)
#  DT                     0.33      2024-04-04 [1] CRAN (R 4.4.0)
#  edgeR                  4.4.1     2024-12-02 [1] Bioconductor 3.20 (R 4.4.2)
#  ExperimentHub          2.14.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  farver                 2.1.2     2024-05-13 [1] CRAN (R 4.4.0)
#  fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.4.0)
#  fields                 16.3      2024-09-30 [1] CRAN (R 4.4.1)
#  filelock               1.0.3     2023-12-11 [1] CRAN (R 4.4.0)
#  forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.4.0)
#  foreach                1.5.2     2022-02-02 [1] CRAN (R 4.4.0)
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.4.0)
#  GenomeInfoDb         * 1.42.1    2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
#  GenomeInfoDbData       1.2.13    2025-01-19 [1] Bioconductor
#  GenomicAlignments      1.42.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  GenomicRanges        * 1.58.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.4.0)
#  ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.4.0)
#  ggrepel              * 0.9.6     2024-09-07 [1] CRAN (R 4.4.1)
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
#  IRanges              * 2.40.1    2024-12-05 [1] Bioconductor 3.20 (R 4.4.2)
#  irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.4.0)
#  iterators              1.0.14    2022-02-05 [1] CRAN (R 4.4.0)
#  jquerylib              0.1.4     2021-04-26 [1] CRAN (R 4.4.0)
#  jsonlite               1.8.9     2024-09-20 [1] CRAN (R 4.4.1)
#  KEGGREST               1.46.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  labeling               0.4.3     2023-08-29 [1] CRAN (R 4.4.0)
#  later                  1.4.1     2024-11-27 [1] CRAN (R 4.4.1)
#  lattice                0.22-6    2024-03-20 [1] CRAN (R 4.4.2)
#  lazyeval               0.2.2     2019-03-15 [1] CRAN (R 4.4.0)
#  lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.4.0)
#  limma                  3.62.2    2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
#  locfit                 1.5-9.10  2024-06-24 [1] CRAN (R 4.4.0)
#  lubridate            * 1.9.4     2024-12-08 [1] CRAN (R 4.4.1)
#  magick                 2.8.5     2024-09-20 [1] CRAN (R 4.4.1)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.4.0)
#  maps                   3.4.2.1   2024-11-10 [1] CRAN (R 4.4.1)
#  Matrix                 1.7-1     2024-10-18 [1] CRAN (R 4.4.2)
#  MatrixGenerics       * 1.18.1    2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
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
#  readxl                 1.4.3     2023-07-06 [1] CRAN (R 4.4.0)
#  rematch2               2.1.2     2020-05-01 [1] CRAN (R 4.4.0)
#  restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.4.0)
#  rjson                  0.2.23    2024-09-16 [1] CRAN (R 4.4.1)
#  rlang                  1.1.5     2025-01-17 [1] CRAN (R 4.4.1)
#  rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.4.0)
#  Rsamtools              2.22.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  RSQLite                2.3.9     2024-12-03 [1] CRAN (R 4.4.1)
#  rstudioapi             0.17.1    2024-10-22 [1] CRAN (R 4.4.1)
#  rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.4.0)
#  rtracklayer            1.66.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  S4Arrays               1.6.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  S4Vectors            * 0.44.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  sass                   0.4.9     2024-03-15 [1] CRAN (R 4.4.0)
#  ScaledMatrix           1.14.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  scales                 1.3.0     2023-11-28 [1] CRAN (R 4.4.0)
#  scater               * 1.34.0    2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
#  scuttle              * 1.16.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.4.0)
#  shiny                  1.10.0    2024-12-14 [1] CRAN (R 4.4.1)
#  shinyWidgets           0.8.7     2024-09-23 [1] CRAN (R 4.4.1)
#  SingleCellExperiment * 1.28.1    2024-11-23 [1] Bioconductor 3.20 (R 4.4.2)
#  spam                   2.11-0    2024-10-03 [1] CRAN (R 4.4.1)
#  SparseArray            1.6.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  SpatialExperiment    * 1.16.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  spatialLIBD          * 1.18.0    2024-11-07 [1] Bioconductor 3.20 (R 4.4.1)
#  statmod                1.5.0     2023-01-06 [1] CRAN (R 4.4.0)
#  stringi                1.8.4     2024-05-06 [1] CRAN (R 4.4.0)
#  stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.4.0)
#  SummarizedExperiment * 1.36.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  systemfonts            1.2.0     2025-01-16 [1] CRAN (R 4.4.1)
#  textshaping            0.4.1     2024-12-06 [1] CRAN (R 4.4.1)
#  tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.4.0)
#  tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.4.0)
#  tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.4.0)
#  tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.4.0)
#  timechange             0.3.0     2024-01-18 [1] CRAN (R 4.4.0)
#  tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.4.0)
#  UCSC.utils             1.2.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  unigd                  0.1.2     2024-06-05 [1] CRAN (R 4.4.0)
#  utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.4.0)
#  vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.4.0)
#  vipor                  0.4.7     2023-12-18 [1] CRAN (R 4.4.0)
#  viridis                0.6.5     2024-01-29 [1] CRAN (R 4.4.0)
#  viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.4.0)
#  vroom                  1.6.5     2023-12-05 [1] CRAN (R 4.4.0)
#  withr                  3.0.2     2024-10-28 [1] CRAN (R 4.4.1)
#  XML                    3.99-0.18 2025-01-01 [1] CRAN (R 4.4.1)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.4.0)
#  XVector                0.46.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
#  yaml                   2.3.10    2024-07-26 [1] CRAN (R 4.4.0)
#  zlibbioc               1.52.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)

#  [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library

# ─────────────────────────────────────────────────────────────────────────────────────────
