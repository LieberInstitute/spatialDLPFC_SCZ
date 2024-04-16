# Load packages ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(scater)
  library(sessioninfo)
})

# Load data ----
## Load spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

## Load outlier keys ----
tot_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

## Load excluded gene df ----
gene_exclude_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "excluded_genes_df.rds"
  )
)

# Prep steps ----
## Remove out_tissue & outlier spots ----
spe$remove <- tot_outlier_df[spe$key, "remove"]
spe <- spe[, spe$remove == FALSE]

## Remove useless colData ----
discard_var <- c(
  "X10x_graphclust", "X10x_kmeans_10_clusters",
  "X10x_kmeans_2_clusters", "X10x_kmeans_3_clusters",
  "X10x_kmeans_4_clusters", "X10x_kmeans_5_clusters",
  "X10x_kmeans_6_clusters", "X10x_kmeans_7_clusters",
  "X10x_kmeans_8_clusters", "X10x_kmeans_9_clusters",
  # QC metrics
  "all_outlier", "sum_umi", "sum_gene", "expr_chrM",
  "expr_chrM_ratio" # ,
  # redundent SPG columns
  # "spg_tissue", "spg_row", "spg_col"
)

for (.var in discard_var) {
  spe[[.var]] <- NULL
}

# Remove 10x sample-wise reducedDim ----
reducedDim(spe, type = "10x_pca") <- NULL
reducedDim(spe, type = "10x_tsne") <- NULL
reducedDim(spe, type = "10x_umap") <- NULL
reducedDimNames(spe)

## Calcualte library size scale factor ----
spe <- scuttle::computeLibraryFactors(spe)
sizeFactors(spe) |> head()

## Remove genes with constant variance ----
kept_gene_index <- !(rownames(spe) %in% rownames(gene_exclude_df))
spe <- spe[kept_gene_index, ]

# which(rowVars(counts(spe)) == 0)
stopifnot(all(rowVars(counts(spe)) != 0))

## Calculate log norm counts ----
spe <- scater::logNormCounts(
  spe,
  size.factors = sizeFactors(spe),
  transform = "log"
)

assayNames(spe)

# Save spe ----
saveRDS(
  spe,
  here(
    "processed-data/rds/02_visium_qc",
    "test_qc_spe_w_spg_N63_no_img.rds"
  )
)


# Session info ----
session_info()
