# SBATCH info
## time 20 min
## mem 10GB

# Load Packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  # library(SpatialExperiment)
  # library(PRECAST)
  library(sessioninfo)
})


## Load PRECAST object -----
PRECASTObj <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_PRECASTObj_semi_inform_k_2-16.rds"
  )
)

# Organize Cluster Results ----
spot_keys <- map(
  PRECASTObj@seulist,
  ~ .x@meta.data$key
)

PRECAST_df <- map(
  .x = PRECASTObj@resList,
  .f = function(res, # RECAST result for each k, as a list
                spot_keys) {
    # Error Prevention: if sample size matches
    if (
      !all(
        sapply(res$cluster, length) == sapply(spot_keys, length)
      )
    ) {
      stop("Unconformatble dimensions")
    }

    ret_vec <- map2(
      .x = res$cluster,
      .y = spot_keys,
      .f = function(vec_cluster, vec_keys) {
        vec_cluster <- sprintf(
          "spd%02d",
          vec_cluster |> as.vector()
        )
        names(vec_cluster) <- vec_keys
        return(vec_cluster)
      }
    ) |> list_c()
    return(ret_vec)
  },
  spot_keys = spot_keys
)

## Organize list to df (n*k) ----
PRECAST_df_final <- do.call(cbind, PRECAST_df)
k_clus <- attr(PRECASTObj@resList, "para_settings")$K # Order of K
colnames(PRECAST_df_final) <- sprintf("PRECAST_%02d", k_clus)


# Save rds ----
saveRDS(
  PRECAST_df_final,
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

# Session info -----
sessioninfo::session_info()

# (Deprecated) Merge with spe object ----
# col_data_df <- PRECAST_df_final2 |>
#   rownames_to_column(var = "key") |>
#   right_join(
#     colData(spe) |> data.frame(),
#     by = c("key"),
#     relationship = "one-to-one"
#   )
# rownames(col_data_df) <- colnames(spe)
# colData(spe) <- DataFrame(col_data_df)
