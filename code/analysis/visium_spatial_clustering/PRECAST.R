suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SpatialExperiment)
  library(Seurat)
  library(PRECAST)
  library(sessioninfo)
})

# Mem requested
# TODO: adjust
# First try 20G


# Config -----------------------------------------------------------
## Path Config -------------------------------------------------------------
fld_data_spatialcluster <- here("processed-data",
                                "rds", "spatial_cluster")
dir.create(
  file.path(fld_data_spatialcluster, "PRECAST"),
  recursive = T, showWarnings = FALSE
)

# TODO: edit the path
file_DLPFC_enrich_csv <- here("code/spatial_clustering/PRECAST",
                              "TableS8_sig_genes_FDR5perc_enrichment.csv")



## Cluster Configuration -------------------------------------------------
n_marker_gene <- 100
# K <- as.numeric(Sys.getenv("SGE_TASK_ID"))
k_min <- 2
k_max <- 16



# Load QC-ed SPE object --------------------------------------------------
spe <- readRDS(here::here(
  "processed-data/rds/",
  #TODO: replace this path
  "test_spe_after_spot_qc_63.rds")
)



# PRECAST Workflow -------------------------------------------------------
## Select spatialDLPFC marker genes  ----------------------------------
gene_df_raw <- read.csv(
  file_DLPFC_enrich_csv
)

gene_df <- gene_df_raw |> 
  filter(spatial_domain_resolution == "Sp09") |> 
  group_by(test) |> 
  arrange(fdr, .by_group = TRUE) |> 
  slice_head(n=n_marker_gene)


stopifnot(all(gene_df$model_type == "enrichment"))
stopifnot(nrow(gene_df) == 9*n_marker_gene)

cat("NOTE (boyiguo1): ",
    gene_df$ensembl |> unique() |> length(),
    " unique gene markers are selected for spatial clustering. \n")

## Create Seurat Object List for PRECAST ----------------------------------
seuList <- unique(spe$sample_id) |> 
  set_names(unique(spe$sample_id)) |> 
  map(.f = function(id) {
    tmp_spe <- spe[, spe$sample_id == id]
    
    tmp_spe$row <- tmp_spe$array_row
    tmp_spe$col <- tmp_spe$array_col
    
    # browser()
    CreateSeuratObject(
      counts=as.matrix(counts(tmp_spe)),
      meta.data=data.frame(colData(tmp_spe)),
      project="PNN")
  })


## Set-up PRECAST -------------------------------------------------------
print("Start PRECAST")
set.seed(1)
PRECASTObj <- CreatePRECASTObject(
  seuList = seuList,
  selectGenesMethod = NULL,
  customGenelist = gene_df$ensembl |> unique()
)

## Add a model setting in advance for a PRECASTObj object.
## verbose =TRUE helps outputing the information in the algorithm.
PRECASTObj <- AddAdjList(
  PRECASTObj,
  platform = "Visium"
)

PRECASTObj <- AddParSetting(
  PRECASTObj, 
  Sigma_equal = FALSE,
  coreNum = ifelse(
    Sys.getenv("SLURM_NTASKS")=="", # on local machine
    8,                              #Boyi's computer. 
    as.numeric(Sys.getenv("SLURM_NTASKS")
    )),
  maxIter = 30,
  verbose = TRUE
)

print("NOTE (boyiguo1): Finish PRECAST set-up")


## Fit PRECAST model -----------------------------------------------------
# tic()
PRECASTObj <- PRECAST(PRECASTObj,
                      K = seq.int(k_min, k_max))
# PRECASTObj <- PRECAST(PRECASTObj,
#                       K = 8)

saveRDS(
  PRECASTObj, 
  file = file.path(
    fld_data_spatialcluster, "PRECAST", 
    paste0("test_PRECASTObj_semi_inform",".rds")
  )
)

print("NOTE (boyiguo1): Finish PRECAST")




## (Deprecated) PRECASTObj integration -----------------------------------
# Necessary step to get cluster from resList
# PRECASTObj <- SelectModel(PRECASTObj)
# seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")


## (Deprecated) Find marker genes ------------------------------------------

# library(Seurat)
# dat_deg <- FindAllMarkers(seuInt)


# (Deprecated) Create final spe object -----------------------------------------
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
# 
# 
# 
# 
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


# Visualize Clustering Result ---------------------------------------------
# # TODO: output
# tmp <- vis_grid_clus(
#   spe,
#   clustervar = "PRECAST_cluster", 
#   return_plots = TRUE
# )


# Session info ------------------------------------------------------------
sessioninfo::session_info()

