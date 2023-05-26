library("here")
library("sessioninfo")

all_metrics <- read.csv(
    here(
        "processed-data",
        "spaceranger",
        "03_merged_metrics",
        "merged_metrics.csv"
    ),
    check.names = FALSE
)
table(all_metrics$study_name)

metrics <-
    subset(
        all_metrics,
        study_name == "spatialDLPFC"
    )
summary(metrics)

#  [cDNA] pg/ul   % Coverage Array Agilent [pg/ul]   Array #          Ave frag length    Brain           cDNA Amp Cycle
# Min.   :  256   Min.   :30.00    Min.   :  716   Length:30          Min.   :360.0   Length:30          Min.   :17.00
# 1st Qu.: 4055   1st Qu.:71.25    1st Qu.: 7690   Class :character   1st Qu.:391.8   Class :character   1st Qu.:17.00
# Median : 7308   Median :80.00    Median :17331   Mode  :character   Median :426.0   Mode  :character   Median :18.00
# Mean   :12867   Mean   :75.17    Mean   :17696                      Mean   :451.0                      Mean   :17.67
# 3rd Qu.:19861   3rd Qu.:80.00    3rd Qu.:24139                      3rd Qu.:481.5                      3rd Qu.:18.00
# Max.   :51825   Max.   :95.00    Max.   :46616                      Max.   :705.0                      Max.   :19.00
# NA's   :4                        NA's   :4                          NA's   :4
# cDNA Input ng           Ct        Est Read Pairs       index_name         index(i7)         index2_workflow_a(i5)
# Min.   :  8.117   Min.   :16.43   Min.   : 75000000   Length:30          Length:30          Length:30
# 1st Qu.: 51.174   1st Qu.:17.07   1st Qu.:178125000   Class :character   Class :character   Class :character
# Median : 94.068   Median :17.50   Median :200000000   Mode  :character   Mode  :character   Mode  :character
# Mean   :138.265   Mean   :17.54   Mean   :197916667
# 3rd Qu.:198.613   3rd Qu.:17.75   3rd Qu.:200000000
# Max.   :518.250   Max.   :19.48   Max.   :318750000
# NA's   :4
# index2_workflow_b(i5)   Sample #          sheet_file          SI cycles      Slide #             Tissue
# Length:30             Length:30          Length:30          Min.   :12.0   Length:30          Length:30
# Class :character      Class :character   Class :character   1st Qu.:14.0   Class :character   Class :character
# Mode  :character      Mode  :character   Mode  :character   Median :14.0   Mode  :character   Mode  :character
#                                                             Mean   :14.3
#                                                             3rd Qu.:15.0
#                                                             Max.   :18.0
#
# Total cDNA ng    slide_serial_capture_area  Sample.ID         Number.of.Spots.Under.Tissue Number.of.Reads
# Min.   : 32.47   Length:30                 Length:30          Min.   :1825                 Min.   : 79655566
# 1st Qu.:166.07   Class :character          Class :character   1st Qu.:3788                 1st Qu.:210222458
# Median :228.15   Mode  :character          Mode  :character   Median :4012                 Median :275278056
# Mean   :306.57                                                Mean   :3960                 Mean   :296318497
# 3rd Qu.:421.59                                                3rd Qu.:4398                 3rd Qu.:310730968
# Max.   :756.80                                                Max.   :4793                 Max.   :700322341
# NA's   :11
# Mean.Reads.per.Spot Mean.Reads.Under.Tissue.per.Spot Fraction.of.Spots.Under.Tissue Median.Genes.per.Spot
# Min.   : 26702      Min.   : 25389                   Min.   :0.3656                 Min.   : 255
# 1st Qu.: 49007      1st Qu.: 45742                   1st Qu.:0.7588                 1st Qu.:1229
# Median : 67211      Median : 63303                   Median :0.8038                 Median :1380
# Mean   : 75065      Mean   : 69824                   Mean   :0.7933                 Mean   :1492
# 3rd Qu.: 78149      3rd Qu.: 73168                   3rd Qu.:0.8811                 3rd Qu.:1697
# Max.   :168347      Max.   :158716                   Max.   :0.9601                 Max.   :2409
#
# Median.UMI.Counts.per.Spot Valid.Barcodes     Valid.UMIs     Sequencing.Saturation Q30.Bases.in.Barcode
# Min.   : 337               Min.   :0.9440   Min.   :0.9961   Min.   :0.8824        Min.   :0.9342
# 1st Qu.:1941               1st Qu.:0.9652   1st Qu.:0.9992   1st Qu.:0.9133        1st Qu.:0.9443
# Median :2257               Median :0.9672   Median :0.9996   Median :0.9417        Median :0.9560
# Mean   :2505               Mean   :0.9659   Mean   :0.9992   Mean   :0.9367        Mean   :0.9523
# 3rd Qu.:3115               3rd Qu.:0.9685   3rd Qu.:0.9998   3rd Qu.:0.9529        3rd Qu.:0.9590
# Max.   :4432               Max.   :0.9721   Max.   :0.9999   Max.   :0.9947        Max.   :0.9669
#
# Q30.Bases.in.RNA.Read Q30.Bases.in.UMI Reads.Mapped.to.Genome Reads.Mapped.Confidently.to.Genome
# Min.   :0.9100        Min.   :0.9370   Min.   :0.7005         Min.   :0.6670
# 1st Qu.:0.9341        1st Qu.:0.9473   1st Qu.:0.8845         1st Qu.:0.8204
# Median :0.9401        Median :0.9541   Median :0.9227         Median :0.8553
# Mean   :0.9381        Mean   :0.9530   Mean   :0.9059         Mean   :0.8429
# 3rd Qu.:0.9478        3rd Qu.:0.9592   3rd Qu.:0.9601         3rd Qu.:0.8918
# Max.   :0.9566        Max.   :0.9657   Max.   :0.9728         Max.   :0.9275
#
# Reads.Mapped.Confidently.to.Intergenic.Regions Reads.Mapped.Confidently.to.Intronic.Regions
# Min.   :0.04495                                Min.   :0.01281
# 1st Qu.:0.05846                                1st Qu.:0.01486
# Median :0.07147                                Median :0.01691
# Mean   :0.07259                                Mean   :0.01812
# 3rd Qu.:0.08195                                3rd Qu.:0.02080
# Max.   :0.10156                                Max.   :0.02964
#
# Reads.Mapped.Confidently.to.Exonic.Regions Reads.Mapped.Confidently.to.Transcriptome Reads.Mapped.Antisense.to.Gene
# Min.   :0.5946                             Min.   :0.5762                            Min.   :0.003124
# 1st Qu.:0.7192                             1st Qu.:0.6975                            1st Qu.:0.004915
# Median :0.7642                             Median :0.7444                            Median :0.006980
# Mean   :0.7522                             Mean   :0.7323                            Mean   :0.007846
# 3rd Qu.:0.7976                             3rd Qu.:0.7801                            3rd Qu.:0.009858
# Max.   :0.8644                             Max.   :0.8374                            Max.   :0.019310
#
# Fraction.Reads.in.Spots.Under.Tissue Total.Genes.Detected summary_file        study_name        study_name_short
# Min.   :0.9421                       Min.   :16946        Length:30          Length:30          Length:30
# 1st Qu.:0.9638                       1st Qu.:21066        Class :character   Class :character   Class :character
# Median :0.9752                       Median :21734        Mode  :character   Mode  :character   Mode  :character
# Mean   :0.9739                       Mean   :21547
# 3rd Qu.:0.9851                       3rd Qu.:22214
# Max.   :0.9925                       Max.   :22841
#
# study_name_sheet   sample_status      tissue_status
# Length:30          Length:30          Length:30
# Class :character   Class :character   Class :character
# Mode  :character   Mode  :character   Mode  :character

## Clean things up for exporting
clean_metrics <- function(info) {
    info$sheet_file <- info$Tissue <- info[["Sample #"]] <- NULL
    info$summary_file <- gsub("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data", "https://github.com/LieberInstitute/spatialDLPFC/tree/main/processed-data", info$summary_file)
    info$study_name <- info$study_name_short <- info$study_name_sheet <- info$sample_status <- info$tissue_status <- NULL
    return(info)
}
head(clean_metrics(metrics))

write.csv(
    clean_metrics(metrics),
    here(
        "processed-data",
        "spaceranger",
        "03_merged_metrics",
        "spatialDLPFC_Visium_metrics.csv"
    )
)

write.csv(
    clean_metrics(subset(
        all_metrics,
        study_name == "spatailDLPFC_IF"
    )),
    here(
        "processed-data",
        "spaceranger",
        "03_merged_metrics",
        "spatialDLPFC_Visium_SPG_metrics.csv"
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 (2022-10-31)
#  os       macOS Ventura 13.0.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/Mexico_City
#  date     2023-02-06
#  rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
#  pandoc   2.17.1.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date (UTC) lib source
#  biocthis      1.9.1   2022-11-01 [1] Github (lcolladotor/biocthis@af38c7c)
#  brio          1.1.3   2021-11-30 [1] CRAN (R 4.2.0)
#  cachem        1.0.6   2021-08-19 [1] CRAN (R 4.2.0)
#  callr         3.7.3   2022-11-02 [1] CRAN (R 4.2.2)
#  cli           3.6.0   2023-01-09 [1] CRAN (R 4.2.0)
#  colorout      1.2-2   2022-03-01 [1] Github (jalvesaq/colorout@79931fd)
#  crayon        1.5.2   2022-09-29 [1] CRAN (R 4.2.0)
#  data.table    1.14.6  2022-11-16 [1] CRAN (R 4.2.0)
#  devtools    * 2.4.5   2022-10-11 [1] CRAN (R 4.2.0)
#  digest        0.6.31  2022-12-11 [1] CRAN (R 4.2.0)
#  dplyr         1.1.0   2023-01-29 [1] CRAN (R 4.2.0)
#  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
#  fansi         1.0.4   2023-01-22 [1] CRAN (R 4.2.0)
#  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
#  fs            1.6.0   2023-01-23 [1] CRAN (R 4.2.0)
#  generics      0.1.3   2022-07-05 [1] CRAN (R 4.2.0)
#  glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
#  here        * 1.0.1   2020-12-13 [1] CRAN (R 4.2.0)
#  hms           1.1.2   2022-08-19 [1] CRAN (R 4.2.0)
#  htmltools     0.5.4   2022-12-07 [1] CRAN (R 4.2.0)
#  htmlwidgets   1.6.1   2023-01-07 [1] CRAN (R 4.2.0)
#  httpuv        1.6.8   2023-01-12 [1] CRAN (R 4.2.0)
#  later         1.3.0   2021-08-18 [1] CRAN (R 4.2.0)
#  lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.2.1)
#  lubridate     1.9.1   2023-01-24 [1] CRAN (R 4.2.0)
#  magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
#  memoise       2.0.1   2021-11-26 [1] CRAN (R 4.2.0)
#  mime          0.12    2021-09-28 [1] CRAN (R 4.2.0)
#  miniUI        0.1.1.1 2018-05-18 [1] CRAN (R 4.2.0)
#  pillar        1.8.1   2022-08-19 [1] CRAN (R 4.2.0)
#  pkgbuild      1.4.0   2022-11-27 [1] CRAN (R 4.2.2)
#  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
#  pkgload       1.3.2   2022-11-16 [1] CRAN (R 4.2.2)
#  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.2.0)
#  processx      3.8.0   2022-10-26 [1] CRAN (R 4.2.0)
#  profvis       0.3.7   2020-11-02 [1] CRAN (R 4.2.0)
#  promises      1.2.0.1 2021-02-11 [1] CRAN (R 4.2.0)
#  prompt        1.0.1   2022-03-01 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps            1.7.2   2022-10-26 [1] CRAN (R 4.2.0)
#  purrr         1.0.1   2023-01-10 [1] CRAN (R 4.2.0)
#  R.cache       0.16.0  2022-07-21 [1] CRAN (R 4.2.0)
#  R.methodsS3   1.8.2   2022-06-13 [1] CRAN (R 4.2.0)
#  R.oo          1.25.0  2022-06-12 [1] CRAN (R 4.2.0)
#  R.utils       2.12.2  2022-11-11 [1] CRAN (R 4.2.0)
#  R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
#  Rcpp          1.0.10  2023-01-22 [1] CRAN (R 4.2.0)
#  remotes       2.4.2   2021-11-30 [1] CRAN (R 4.2.0)
#  rlang         1.0.6   2022-09-24 [1] CRAN (R 4.2.0)
#  rprojroot     2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
#  rsthemes      0.3.1   2022-03-01 [1] Github (gadenbuie/rsthemes@bbe73ca)
#  rstudioapi    0.14    2022-08-22 [1] CRAN (R 4.2.0)
#  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
#  shiny         1.7.4   2022-12-15 [1] CRAN (R 4.2.2)
#  stringi       1.7.12  2023-01-11 [1] CRAN (R 4.2.0)
#  stringr       1.5.0   2022-12-02 [1] CRAN (R 4.2.0)
#  styler        1.9.0   2023-01-15 [1] CRAN (R 4.2.0)
#  suncalc       0.5.1   2022-09-29 [1] CRAN (R 4.2.0)
#  testthat    * 3.1.6   2022-12-09 [1] CRAN (R 4.2.0)
#  tibble        3.1.8   2022-07-22 [1] CRAN (R 4.2.1)
#  tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.2.0)
#  timechange    0.2.0   2023-01-11 [1] CRAN (R 4.2.0)
#  urlchecker    1.0.1   2021-11-30 [1] CRAN (R 4.2.0)
#  usethis     * 2.1.6   2022-05-25 [1] CRAN (R 4.2.0)
#  utf8          1.2.3   2023-01-31 [1] CRAN (R 4.2.0)
#  vctrs         0.5.2   2023-01-23 [1] CRAN (R 4.2.0)
#  xtable        1.8-4   2019-04-21 [1] CRAN (R 4.2.0)
#
#  [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
