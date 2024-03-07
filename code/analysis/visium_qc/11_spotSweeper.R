# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SpotSweeper)
  library(sessioninfo)
  library(escheR)
  library(here)
  library(scater)
  # library(spatialLIBD)
  library(tidyverse)
})

# File Paths --------------------------------------------------------------
path_raw_spe <- here(
  "processed-data/rds",
  "01_build_spe",
  # "test_spe_w_spg_n4.rds" # Test only
  "raw_spe_wo_SPG_N63.rds"
)

# Load SPE -----------------------------------------------------------------
raw_spe <- readRDS(
  path_raw_spe
)

# in-tissue spots only
raw_spe <- raw_spe[, raw_spe$in_tissue == TRUE]

# Create logcounts
if (!"logcounts" %in% assayNames(raw_spe)) {
  raw_spe <- raw_spe[, raw_spe$sum_umi != 0]
  raw_spe <- logNormCounts(raw_spe)
}





# Spot Sweeper

## Find Tissue Artifacts (Tears and Hangneils)---------------------------------------------------

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

# pdf(
#   here("plots/02_visium_qc/spotSweeper", "artifacts_tears.pdf"),
#   width = 5.9, height = 3.9
# )
for (.sample in artf_samples) {
  # .sample <- artf_samples[[1]]
  sub_spe <- raw_spe[, raw_spe$sample_id == .sample]
  sub_spe <- sub_spe[, sub_spe$in_tissue == 1]

  sub_spe <- findArtifacts(
    sub_spe,
    mito_percent = "expr_chrM_ratio",
    mito_sum = "expr_chrM",
    n_rings = 5,
    name = "artifact"
  )

  # print(
  #   make_escheR(sub_spe) |>
  #     add_fill(var = "expr_chrM_ratio", point_size = 1.25) |>
  #     add_ground(var = "artifact", point_size = 1.25, stroke = 0.3) +
  #     scale_fill_gradient(low = "white", high = "black") +
  #     scale_color_manual(
  #       name = "Affected area", # turn off legend name for ground_truth
  #       values = c(
  #         "TRUE" = "red",
  #         "FALSE" = "transparent"
  #       )
  #     ) +
  #     labs(title = .sample)
  # )
  # Warning message:
  #   In .check_reddim_names(x, value, withDimnames) :
  #   non-NULL 'rownames(value)' should be the same as 'colnames(x)' for
  # 'reducedDim<-'. This will be an error in the next release of Bioconductor.

  artf_spots_key <- c(artf_spots_key, sub_spe$key[sub_spe$artifact == TRUE])
}
# dev.off()

raw_spe$scratch <- "No"
raw_spe$scratch[raw_spe$key %in% artf_spots_key] <- "Yes"
saveRDS(
  artf_spots_key,
  here("processed-data/visium_qc/scratch_key.rds"),
)

# Spot-level ----
# raw_spe <- raw_spe[, !(raw_spe$key %in% artf_spots_key)]
# spe_clean_artifacts <- spe_clean_artifacts[, spe_clean_artifacts$in_tissue == TRUE]


# pdf(
#   here("plots/02_visium_qc/spotSweeper", "outliers_n_ngb36.pdf"),
#   width = 5.9, height = 3.9
# )
outlier_spots_key <- c()
for (.sample in unique(raw_spe$sample_id)) {
  # .sample <- unique(raw_spe$sample_id)[1]

  sub_spe <- raw_spe[, raw_spe$sample_id == .sample]


  features <- c("sum_umi", "sum_gene")
  sub_spe <- localOutliers(
    sub_spe,
    features = features,
    # n_neighbors=18,
    n_neighbors = 36,
    data_output = TRUE,
    method = "multivariate"
  )

  outlier_spots_key <- c(
    outlier_spots_key,
    sub_spe$key[sub_spe$local_outliers == TRUE]
  )
  #   print(
  #     make_escheR(sub_spe) |>
  #       add_fill(var = "sum_umi_log2", point_size = 1.25) |>
  #       add_ground(var = "local_outliers", point_size = 1.25, stroke = 1) +
  #       scale_color_manual(
  #         name = "", # turn off legend name for ground_truth
  #         values = c(
  #           "TRUE" = "red",
  #           "FALSE" = "transparent"
  #         )
  #       ) +
  #       scale_fill_gradient(low = "white", high = "darkgreen")
  #   )


  #   #
  #   # make_escheR(sub_spe) |>
  #   #   add_fill(var = "sum_gene_log2", point_size=1.25) |>
  #   #   add_ground(var = "local_outliers",point_size=1.25, stroke = 1) +
  #   #   scale_color_manual(
  #   #     name = "", # turn off legend name for ground_truth
  #   #     values = c(
  #   #       "TRUE" = "red",
  #   #       "FALSE" = "transparent")
  #   #   ) +
  #   #   scale_fill_gradient(low ="white",high =  "darkgreen")
  #   #
  #   # make_escheR(sub_spe) |>
  #   #   add_fill(var = "sum_umi_z", point_size=1.25) +
  #   #   # add_ground(var = "local_outliers",point_size=1.25, stroke = 1) +
  #   #   # scale_color_manual(
  #   #   #   name = "", # turn off legend name for ground_truth
  #   #   #   values = c(
  #   #   #     "TRUE" = "red",
  #   #   #     "FALSE" = "transparent")
  #   #   # ) +
  #   #   scale_fill_gradient(low ="white",high =  "darkgreen")
}
# dev.off()


raw_spe$local_outliers <- "FALSE"
raw_spe$local_outliers[raw_spe$key %in% outlier_spots_key] <- "TRUE"

saveRDS(
  outlier_spots_key,
  here("processed-data/visium_qc/outlier_key_ngb_36.rds")
)

# Validation Plots -----
## Dimension Reduction (PCA & UMAP) -------------------------------------------------------------
set.seed(1)
raw_spe <- BiocSingular::runPCA(raw_spe) # 50 PCs
raw_spe <- scater::runUMAP(raw_spe, dimred = "PCA")

# Note:
#     * t-SNE is computationally infeasible for the large dataset
#     * UMAP is merely as visual check.

## Plots
# plotReducedDim(
#   spe_dimRed,
#   dimred = "UMAP",
#   color_by = "sample_id")

## Mito_ratio vs sum_gene (miQC) validation ----
max_mito_ratio <- max(raw_spe$expr_chrM_ratio)
max_sum_gene <- max(raw_spe$sum_gene)

## Vis mito_ratio vs sum_gene relationship
pdf(here("plots/02_visium_qc/spotSweeper/test_validation_mt_gene_relationship_scratch.pdf"))
# for (.sample in unique(raw_spe$sample_id)) {
for (.sample in artf_samples) {
  # .sample <- unique(raw_spe$sample_id)[1]

  # TODO: why there's warning
  # Warning messages:
  # 1: Removed 31 rows containing missing values (`geom_point()`).
  # 2: Removed 13 rows containing missing values (`geom_point()`).
  ret_p <- ggplot() +
    geom_point(
      aes(
        x = raw_spe[, raw_spe$sample_id == .sample]$sum_gene,
        y = raw_spe[, raw_spe$sample_id == .sample]$expr_chrM_ratio,
        color = raw_spe[, raw_spe$sample_id == .sample]$scratch
      ),
      size = 0.3
    ) +
    geom_rug(
      aes(
        x = raw_spe[, raw_spe$sample_id == .sample]$sum_gene,
        y = raw_spe[, raw_spe$sample_id == .sample]$expr_chrM_ratio,
        color = raw_spe[, raw_spe$sample_id == .sample]$scratch
      ),
      size = 0.3
    ) +
    labs(title = .sample) +
    scale_color_discrete(name = "scratch") +
    scale_x_continuous(limits = c(0, max_sum_gene)) +
    scale_y_continuous(limits = c(0, max_mito_ratio))

  print(ret_p)
}
dev.off()


pdf(here("plots/02_visium_qc/spotSweeper/test_validation_mt_gene_relationship_outlier.pdf"))
# for (.sample in unique(raw_spe$sample_id)) {
for (.sample in unique(raw_spe$sample_id)) {
  # .sample <- unique(raw_spe$sample_id)[1]

  # TODO: why there's warning
  # Warning messages:
  ret_p <- ggplot() +
    geom_point(
      aes(
        x = raw_spe[, raw_spe$sample_id == .sample]$sum_gene,
        y = raw_spe[, raw_spe$sample_id == .sample]$expr_chrM_ratio,
        color = raw_spe[, raw_spe$sample_id == .sample]$local_outliers
      ),
      size = 0.3
    ) +
    geom_rug(
      aes(
        x = raw_spe[, raw_spe$sample_id == .sample]$sum_gene,
        y = raw_spe[, raw_spe$sample_id == .sample]$expr_chrM_ratio,
        color = raw_spe[, raw_spe$sample_id == .sample]$local_outliers
      ),
      size = 0.3
    ) +
    labs(title = .sample) +
    scale_color_discrete(name = "local_outliers") +
    scale_x_continuous(limits = c(0, max_sum_gene)) +
    scale_y_continuous(limits = c(0, max_mito_ratio))

  print(ret_p)
}
dev.off()


# Session Info ------------------------------------------------------------
sessionInfo()
