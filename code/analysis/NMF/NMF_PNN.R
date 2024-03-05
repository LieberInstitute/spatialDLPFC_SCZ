# Load Library ------------------------------------------------------------
suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
    library(sessioninfo)
    library(tidyverse)
    library(spatialLIBD)
    library(tidyverse)
})


# Load Data -----
## Load NMF Res --------
model <- readRDS(
    here(
        "processed-data/rds/NMF/",
        "test_NMF_all_k50.rds"
    )
)

## Load Spe ---------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

spe <- readRDS(
  path_PRECAST_int_spe
)

spot_key <- spe$key



raw_spe <- readRDS(
    here::here(
        "processed-data", "rds", "01_build_spe",
        # TODO: rename
        "test_raw_spe_w_spg_N63.rds"
    )
)

spe <- raw_spe[, raw_spe$key %in% spot_key]



# For test only
# spe_backup <- spe
# spe <- spe_backup

spe$dx <- metadata(spe)$dx_df$dx[
    match(
        spe$sample_id,
        metadata(spe)$dx_df$sample_id
    )
]

# k <- 2
k <- 50
patterns <- t(model$h) # these are the factors
colnames(patterns) <- paste0("NMF", 1:k)
colData(spe) <- cbind(colData(spe), patterns)

spe_complete <- spe[, !is.na(spe$spg_PBW)]



# Overall NMF Correlation ----
cor(
    colData(spe_complete) |> data.frame() |> select(starts_with("NMF")),
    spe_complete$spg_PBW
)

sample_nmf_mat <- list()
# NMF_PNN Corr per sample ----
for (.sample in unique(spe_complete$sample_id)) {
    # .sample <- unique(spe_complete$sample_id)[1]
    ret <- cor(
        colData(spe_complete[, spe_complete$sample_id == .sample]) |>
            data.frame() |> select(starts_with("NMF")),
        spe_complete[, spe_complete$sample_id == .sample]$spg_PBW
    )
    sample_nmf_mat[[.sample]] <- ret 
}

sample_nmf_mat <- do.call(cbind, sample_nmf_mat)
colnames(sample_nmf_mat) <- unique(spe_complete$sample_id)



heatmap(sample_nmf_mat, scale = "none")

heatmap(abs(sample_nmf_mat), scale = "none")


# Interpret NMF 5
str(model$w)
NMF5_loadings <- model$w[,5]
names(NMF5_loadings) <- rowData(spe)$gene_name


NMF5_loadings |> sort( decreasing = TRUE ) |> head(100) |> names()

which((NMF5_loadings |> sort( decreasing = TRUE ) |> names() )== "PVALB")

# Session Info ----
sessioninfo::session_info()