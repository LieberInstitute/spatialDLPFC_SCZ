# SBATCH info
## time 20 min
## mem 30GB

# Load Packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  # library(SpatialExperiment)
  library(PRECAST)
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
PRECAST_df_final <- PRECAST_df_final |>
  data.frame() |>
  rownames_to_column(var = "key")




# Save Data ----
## PRECAST labels as data.frame ----
saveRDS(
  PRECAST_df_final,
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
  # NOTE: later manually renamed to
  # here(
  #   "processed-data/rds/03_visium_spatial_clustering",
  #   "PRECAST_label_df_semi_sup_k_2-16.rds"
  # )
)


## spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)

col_data_df <- PRECAST_df_final |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

saveRDS(
  spe,
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# Deprecated code ----
##  Load PRECAST data ----
# PRECASTObj_model_select <- SelectModel(PRECASTObj)

# # NOTE:
# # Boyi deoesn't really like this step, in the sense that the batch correction step is slightly iffy for the following reasons
# # 1. threatically speaking, the batch correction is done with house keeping genes.
# # How often the author would update the house keeping genes and where to acquire it.

# # Boyi's Recommentation
# # It would probably better to export the PRECAST latent embedding, and run batch correction methods external to this step.

# # Note:  Integration based on PRECAST pipeline (not used.)
# seuInt <- IntegrateSpaData(PRECASTObj_model_select, species = "Human")

# seuInt <- AddUMAP(seuInt,
#   seed = 1
# )
# seuInt <- AddTSNE(seuInt, seed = 1)

# seuInt <- readRDS(
#   here(
#     "processed-data/rds/spatial_cluster/PRECAST",
#     "test_seuInt_UMAP_tsne.rds"
#   )
# )

# Session info -----
sessioninfo::session_info()
