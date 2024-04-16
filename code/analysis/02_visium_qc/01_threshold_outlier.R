# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(scater)
  library(tidyverse)
  library(SpotSweeper)
  library(sessioninfo)
  library(here)
})

stopifnot(
  "To replicate the result, please use spotsweeper 0.99.2" =
    package.version("SpotSweeper") >= "0.99.2"
)

# Load Data ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

## Remove out tissue spots ----
ret_spe <- spe[, spe$in_tissue == TRUE]

# Spot Annotation ----
## Annotate compromised regions (SpotSweeper)----
sample_info <- readxl::read_xlsx(
  here(
    "raw-data/experiment_info",
    "PNN_64_collection_Summary_final.xlsx"
  )
)

artf_samples <- sample_info |>
  select(SlideSerial, CaptureArea, hangnail, tear) |>
  filter(tear == "yes") |>
  transmute(sample_id = paste0(SlideSerial, "_", CaptureArea)) |>
  pull(sample_id)


artf_spots_key <- c()

for (.sample in artf_samples) {
  sub_spe <- ret_spe[, ret_spe$sample_id == .sample]
  sub_spe <- sub_spe[, sub_spe$in_tissue == TRUE]

  # Prevent subset mis-behaving
  stopifnot(ncol(sub_spe) != 0)

  sub_spe <- findArtifacts(
    sub_spe,
    mito_percent = "expr_chrM_ratio",
    mito_sum = "expr_chrM",
    n_rings = 5,
    name = "artifact"
  )

  artf_spots_key <- c(artf_spots_key, sub_spe$key[sub_spe$artifact == TRUE])
}


# Summarize outlier df ----
outlier_df <- data.frame(
  key = spe$key,
  out_tissue_spots = !spe$in_tissue,
  umi_lt_100 = spe$sum_umi <= 100,
  gene_lt_200 = spe$sum_gene <= 200,
  artifact = (spe$key %in% artf_spots_key)
)

saveRDS(
  outlier_df,
  file = here::here(
    "processed-data/rds/02_visium_qc",
    "outlier_df.rds"
  )
)

# Session info ----
session_info()
