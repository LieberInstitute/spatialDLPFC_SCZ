# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SpatialExperiment)
  library(PRECAST)
  library(sessioninfo)
})


# Path --------------------------------------------------------------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_precast_output <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_PRECASTObj_semi_inform",".rds")
)

path_spe <- file.path(
  here::here(
    "processed-data/rds/",
    #TODO: replace this path
    "test_spe_after_spot_qc_63.rds")
)

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)


# Load PRECAST Data -------------------------------------------------------
PRECASTObj <- readRDS(
  path_precast_output
)

# Organize Cluster Results ------------------------------------------------

k_clus <- attr(PRECASTObj@resList, "para_settings")$K

# Assign name to each resList
precast_resList <- map2(
  PRECASTObj@seulist,
  PRECASTObj@resList,
  .f = function(seu_obj, res_obj){
    browser()
    names(seu_obj)
  }
)

sample_names <- names(PRECASTObj@seulist)
sample_sizes <- sapply(PRECASTObj@seulist, FUN = ncol)

spot_keys <- map(PRECASTObj@seulist, ~.x@meta.data$key)


PRECAST_df <- map2(
  .x = PRECASTObj@resList,
  .y = k_clus,
  .f = function(
    res,      # RECAST result for each k, as a list
    k,         # Number of Clusters,
    spot_keys
  ){

    # Error Prevention: if sample size matches
    if(
      !all(
        sapply(res$cluster, length) == sapply(spot_keys, length)
      )
    ){
      stop("Unconformatble dimensions")
    }
      
    # browser()

  ret_vec <- map2(
    .x = res$cluster,
    .y = spot_keys,
    .f = function(vec_cluster, vec_keys){
      # browser()
      vec_cluster <- vec_cluster |> as.vector() |> factor()
      names(vec_cluster) <- vec_keys
      return(vec_cluster)
    }
  ) |> list_c()
  # browser()
  # ret_vec <- paste("Clus", ret_vec) |> factor()
  
    ret_vec

  },
  spot_keys = spot_keys
)

# Save to spe -------------------------------------------------------------
# spe <- readRDS(
# )

PRECAST_df_final <- do.call(cbind, PRECAST_df)
colnames(PRECAST_df_final) <- paste0("PRECAST_", k_clus)
PRECAST_df_final2 <- PRECAST_df_final |> as.data.frame()


spe <- readRDS(path_spe)

col_data_df <- PRECAST_df_final2 |>
   rownames_to_column(var = "key") |> 
   right_join(
     colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)


saveRDS(
  spe,
  file =
    path_PRECAST_int_spe
)


# SessionInfo -------------------------------------------------------------
sessioninfo::session_info()

