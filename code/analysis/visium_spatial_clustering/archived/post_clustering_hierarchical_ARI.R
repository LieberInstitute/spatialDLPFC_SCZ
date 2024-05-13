# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  # library(spatialLIBD)
  library(tidyverse)
  # library(escheR)
  # library(ggpubr)
  library(mclust)
  library(sessioninfo)
  
})

# Path --------------------------------------------------------------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

# Load data ---------------------------------------------------------------

spe <- readRDS(
  path_PRECAST_int_spe
)

col_df <- colData(spe) |> data.frame()

# spe_backup <- spe


# Pairwise_ARI ------------------------------------------------------------
# .sample_id <- spe$sample_id[1]

ARI_list <- list()

for(.sample_id in unique(spe$sample_id)){
  spe_sub <- spe[, spe$sample_id == .sample_id]
  sample_ari <- seq.int(2, 11) |> 
    map_dbl(.f = function(k){
      adjustedRandIndex(
        spe_sub[[paste0("PRECAST_", k)]],
        spe_sub[[paste0("PRECAST_", k+1)]]
        
      )
    })
  
  ARI_list[[.sample_id]] <- sample_ari |> t() |> data.frame()
} 


ARI_mat <- ARI_list |> list_rbind()

rownames(ARI_mat) <- names(ARI_list)

ARI_mat[,1] |> sort()

# ARI_mat <- ARI_list |> list_c() |> 
#   matrix(ncol = seq.int(2, 11) |> length(), byrow = TRUE)
# 
# colnames(ARI_mat) <- paste0("k_", seq.int(2, 11))
# rownames(ARI_mat) <- unique(spe$sample_id)

heatmap(
  ARI_mat, Colv = NA, Rowv = NA, 
  # scale = "none"
  scale = "column"
)


# High ARI: didn't add much space in the cluster
# Low ARI: a lot of areas either being new cluster (with substaintial spatial area) or cluster differently.

# tmp_backup <- ARI_list
# 




