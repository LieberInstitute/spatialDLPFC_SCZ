# Load packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SpatialExperiment)
  library(Seurat)
  library(PRECAST)
  library(sessioninfo)
})

# Load data
## Load spe object ----
spe <- readRDS(
  here::here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)

# Prep supervised marker genes ----
## Select spatialDLPFC marker genes  -----
gene_df_raw <- read.csv(
  here(
    # TODO: change path
    "code/analysis/visium_spatial_clustering",
    "TableS8_sig_genes_FDR5perc_enrichment.csv"
  )
)

n_marker_gene <- 100
gene_df <- gene_df_raw |>
  filter(spatial_domain_resolution == "Sp09") |>
  group_by(test) |>
  arrange(fdr, .by_group = TRUE) |>
  slice_head(n = n_marker_gene)

stopifnot(all(gene_df$model_type == "enrichment"))
stopifnot(nrow(gene_df) == 9 * n_marker_gene)

cat(
  "NOTE (boyiguo1): ",
  gene_df$ensembl |> unique() |> length(),
  " unique gene markers are selected for spatial clustering. \n"
)

write_csv(
  gene_df,
  here(
    "manuscript/tables",
    "stable_selected_gene_for_PRECAST.csv"
  )
)



# PRECAST Workflow -----
## Create Seurat List for PRECAST -----
seuList <- unique(spe$sample_id) |>
  # seuList <- unique(spe$sample_id)[1:4] |> # Testing
  set_names(unique(spe$sample_id)) |>
  map(.f = function(id) {
    tmp_spe <- spe[, spe$sample_id == id]

    tmp_spe$row <- tmp_spe$array_row
    tmp_spe$col <- tmp_spe$array_col

    CreateSeuratObject(
      counts = as.matrix(counts(tmp_spe)),
      meta.data = data.frame(colData(tmp_spe)),
      project = "PNN"
    )
  })


## Set-up PRECAST -----
print("Start PRECAST")
set.seed(20240416)
PRECASTObj <- CreatePRECASTObject(
  seuList = seuList,
  selectGenesMethod = NULL,
  customGenelist = gene_df$ensembl |> unique(),
  # Set to 0s to force not removing anything
  premin.spots = 0,
  premin.features = 0,
  postmin.spots = 0,
  postmin.features = 0,
  rawData.preserve = FALSE,
  verbose = TRUE
)

# Save mem
rm(spe)
rm(seuList)
gc(verbose = FALSE)

## Add a model setting in advance for a PRECASTObj object.
## verbose =TRUE helps outputing the information in the algorithm.
PRECASTObj <- AddAdjList(
  PRECASTObj,
  type = "fixed_distance",
  platform = "Visium"
)

PRECASTObj <- AddParSetting(
  PRECASTObj,
  Sigma_equal = FALSE,
  coreNum = ifelse(
    Sys.getenv("SLURM_NTASKS") == "", # on local machine
    8, # Boyi's computer.
    as.numeric(Sys.getenv("SLURM_NTASKS"))
  ),
  maxIter = 30,
  verbose = TRUE,
  seed = 20240416
)

print("NOTE (boyiguo1): Finish PRECAST set-up")


## Fit PRECAST model ----
k_min <- 2
k_max <- 16
PRECASTObj <- PRECAST(
  PRECASTObj,
  K = seq.int(k_min, k_max),
  q = 15 # Arbitary number
)

# Save PRECAST Obj ----
saveRDS(
  PRECASTObj,
  file =
    here(
      "processed-data/rds/spatial_cluster",
      "PRECAST",
      "test_PRECASTObj_semi_inform_k_2-16.rds"
    )
)

print("NOTE (boyiguo1): Finish PRECAST")


# Session info ----
sessioninfo::session_info()

# (Deprecated) Code -----
## (Deprecated) PRECASTObj integration -----------------------------------
# Necessary step to get cluster from resList
# PRECASTObj <- SelectModel(PRECASTObj)
# seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
## (Deprecated) Find marker genes ------------------------------------------
# library(Seurat)
# dat_deg <- FindAllMarkers(seuInt)
## (Deprecated) Create final spe object -----------------------------------------
# PRECASTObj <- readRDS(
#   file.path(
#     fld_data_spatialcluster, "PRECAST",
#     paste0("test_PRECASTObj_semi_inform",".rds")
#   )
# )
#
# # Number of cluster
# attr(PRECASTObj@resList, "para_settings")$K
#
# # TODO: write a wrapper function to read in results
# PRECASTObj@resList[[1]]$cluster
# col_data_df <- seuInt@meta.data |>
#   mutate(cluster = factor(cluster)) |>
#   rename_with(~ paste0("PRECAST_", .x)) |>
#   rownames_to_column(var = "key") |>
#   right_join(
#     colData(spe) |> data.frame(),
#     by = c("key"),
#     relationship = "one-to-one"
#   )
# rownames(col_data_df) <- colnames(spe)
# colData(spe) <- DataFrame(col_data_df)
#
# fld_data_spatialcluster <- here("processed-data", "rds", "spatial_cluster")
# dir.create(
#   file.path(fld_data_spatialcluster, "PRECAST"),
#   recursive = T, showWarnings = FALSE
# )
# saveRDS(
#   seuInt,
#   file = file.path(
#     fld_data_spatialcluster, "PRECAST",
#     paste0("test_seuIntObj_semi_inform_K",K,".rds")
#   )
# )
# spe <- readRDS(here(
#   file = file.path(
#     fld_data_spatialcluster, "PRECAST",
#     paste0("test_PRECASTObj_semi_inform_K",K,".rds")
#   )))
## Visualize Clustering Result ----
# # TODO: output
# tmp <- vis_grid_clus(
#   spe,
#   clustervar = "PRECAST_cluster",
#   return_plots = TRUE
# )
