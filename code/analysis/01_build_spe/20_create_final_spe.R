# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
})


# Load data ----
## load raw data ----
spe <- readRDS(
  here::here(
    "processed-data/rds/01_build_spe",
    "test_raw_spe_w_spg_N63.rds"
    # Or using the with image version
    # "raw_spe_w_SPG_N63_loaded_img.rds"
  )
)


# Manipulate colData ----
## Add donor informaiton ----
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]
spe$DX <- toupper(spe$dx)

spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$age <- metadata(spe)$dx_df$age[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$rin <- metadata(spe)$dx_df$RIN[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$slide_id <- sapply(strsplit(spe$sample_id, "_"), function(x) x[1])
spe$lot_num <- sapply(strsplit(spe$slide_id, "-"), function(x) x[1])

spe$brnum <- metadata(spe)$dx_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sample_label <- spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

## Add outlier information -----
tot_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

spe$outlier_cat <- tot_outlier_df[spe$key, "all_outlier"]
# Categories:
# local: spotsweeper local outlier
# outlier: regional artifact + (umi < 100 or gene <= 200)
# Both: marked as outlier using both approach local neither outlier
# neither: not outlier

spe$removed_spots <- tot_outlier_df[spe$key, "remove"]
# NOTE: removed spots include outlier and non-tissue spots
table(spe$removed_spots)

## Add spatial domain information ----
PRECAST_df <- readRDS(
  here(
    "processed-data/rds/03_visium_spatial_clustering",
    "PRECAST_label_df_semi_sup_k_2-16.rds"
  )
)

### Merge PRECAST df to colData ----
# precast_vars <- grep(
#   "^PRECAST_", colnames(PRECAST_df),
#   value = TRUE
# )

col_data_df <- PRECAST_df |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

table(col_data_df$PRECAST_07, useNA = "always")
# Match with outlier N_outlier = 34,689

rownames(col_data_df) <- col_data_df$key

# error prevention
stopifnot(nrow(col_data_df) == ncol(spe))

# DataFrame(col_data_df)|> names()

# NOTE: not sure why this line doesn't work any more
# colData(spe) <- col_data_df |> select(-sample_id) |> DataFrame()
# rownames(colData(spe)) <- spe$key

precast_vars <- grep("^PRECAST", colnames(PRECAST_df), value = TRUE)

# Creates mismatch in merged results
# for (var in precast_vars) {
#   spe[[var]][spe$key] <- col_data_df[, var]
# }

rownames(col_data_df)


# spe$PRECAST_07 <- col_data_df$PRECAST_07
spe$PRECAST_02 <- col_data_df$PRECAST_02[match(spe$key, col_data_df$key)]
spe$PRECAST_03 <- col_data_df$PRECAST_03[match(spe$key, col_data_df$key)]
spe$PRECAST_04 <- col_data_df$PRECAST_04[match(spe$key, col_data_df$key)]
spe$PRECAST_05 <- col_data_df$PRECAST_05[match(spe$key, col_data_df$key)]
spe$PRECAST_06 <- col_data_df$PRECAST_06[match(spe$key, col_data_df$key)]
spe$PRECAST_07 <- col_data_df$PRECAST_07[match(spe$key, col_data_df$key)]
spe$PRECAST_08 <- col_data_df$PRECAST_08[match(spe$key, col_data_df$key)]
spe$PRECAST_09 <- col_data_df$PRECAST_09[match(spe$key, col_data_df$key)]
spe$PRECAST_10 <- col_data_df$PRECAST_10[match(spe$key, col_data_df$key)]
spe$PRECAST_11 <- col_data_df$PRECAST_11[match(spe$key, col_data_df$key)]
spe$PRECAST_12 <- col_data_df$PRECAST_12[match(spe$key, col_data_df$key)]
spe$PRECAST_13 <- col_data_df$PRECAST_13[match(spe$key, col_data_df$key)]
spe$PRECAST_14 <- col_data_df$PRECAST_14[match(spe$key, col_data_df$key)]
spe$PRECAST_15 <- col_data_df$PRECAST_15[match(spe$key, col_data_df$key)]
spe$PRECAST_16 <- col_data_df$PRECAST_16[match(spe$key, col_data_df$key)]

# error prevention
stopifnot(all(precast_vars %in% names(colData(spe))))

stopifnot(
  identical(
    table(PRECAST_df$PRECAST_07),
    table(spe$PRECAST_07)
  )
)

stopifnot(
  identical(
    table(spe$sample_id, spe$PRECAST_07),
    table(col_data_df$sample_id, col_data_df$PRECAST_07)
  )
)



### create annotated spd labels (PRECAST_07) ----
spe$fnl_spd <- spe$PRECAST_07
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(
    anno_lab = factor(
      paste0(gsub("spd", "SpD", spd), "-", label),
      levels = c(
        "SpD07-L1",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      )
    )
  )

spe$spd_label <- spd_anno_df[
  match(spe$fnl_spd, spd_anno_df$spd),
  "anno_lab"
]

# error prevention
# Note: the numbers should be the same as the PRECAST_07 table
table(spe$spd_label, useNA = "always")


# Log transformation ----
## Identify spots that can't be log normalized -----
spe <- scuttle::computeLibraryFactors(spe)

# Note 18 spots doesn't have any RNA sequenced
sum(sizeFactors(spe) == 0)
# [1] 18

non_zero_keys <- which(sizeFactors(spe) != 0) |> names()

spe <- spe[, spe$key %in% non_zero_keys]

## Log normalized -----
spe <- scater::logNormCounts(
  spe,
  size.factors = sizeFactors(spe),
  transform = "log"
)


# Currently deprecated ----
## Remove useless colData ----
# discard_var <- c(
#   # 10x output
#   "X10x_graphclust", "X10x_kmeans_10_clusters",
#   "X10x_kmeans_2_clusters", "X10x_kmeans_3_clusters",
#   "X10x_kmeans_4_clusters", "X10x_kmeans_5_clusters",
#   "X10x_kmeans_6_clusters", "X10x_kmeans_7_clusters",
#   "X10x_kmeans_8_clusters", "X10x_kmeans_9_clusters",
#   "10x_graphclust", "10x_kmeans_10_clusters",
#   "10x_kmeans_2_clusters", "10x_kmeans_3_clusters",
#   "10x_kmeans_4_clusters", "10x_kmeans_5_clusters",
#   "10x_kmeans_6_clusters", "10x_kmeans_7_clusters",
#   "10x_kmeans_8_clusters", "10x_kmeans_9_clusters",
#   # QC metrics
#   "all_outlier", "sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio", "remove" # ,
#   # redundent SPG columns
#   # "spg_tissue", "spg_row", "spg_col"
# )

# for (.var in discard_var) {
#   spe[[.var]] <- NULL
# }

## Remove 10x sample-wise reducedDim ----
# reducedDim(spe, type = "10x_pca") <- NULL
# reducedDim(spe, type = "10x_tsne") <- NULL
# reducedDim(spe, type = "10x_umap") <- NULL
# reducedDimNames(spe)


# Create RDS file ----
saveRDS(
  spe,
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_all_spots.rds"
  )
)

## Create outlier free spe ----
saveRDS(
  spe[, spe$removed_spots == FALSE],
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

# Session info ----
sessioninfo::session_info()
