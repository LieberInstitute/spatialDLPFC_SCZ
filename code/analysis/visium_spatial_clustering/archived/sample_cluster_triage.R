library(here)
library("dplyr")
library("purrr")
library("Seurat")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("PRECAST")
# library(glue)


# Mem requested
# First try 20G


# Path Config -------------------------------------------------------------

fld_data_spatialcluster <- here("processed-data", "rds", "spatial_cluster")
dir.create(
  file.path(fld_data_spatialcluster, "PRECAST"),
  recursive = T, showWarnings = FALSE
)



# Load QC-ed SPE object --------------------------------------------------
raw_spe <- readRDS(here::here(
  "processed-data/rds/",
  #TODO: replace this path
  "test_spe_after_spot_qc_63.rds")
)


bad_sample_name_df <- read.csv(
  here("code/analysis/visium_spatial_clustering",
       "visiumspg_pnn_triage_spatialqc.csv")
) |> filter(VisualCheckSR != "ok")

nrow(bad_sample_name_df)

bad_sampe_name <- with(bad_sample_name_df,
                       paste( SlideSerial, Array,sep="_"))




spe <- raw_spe[, raw_spe$sample_id %in% bad_sampe_name]
stopifnot(
  identical(
    unique(spe$sample_id),
    bad_sampe_name)
)

# Run with all samples
# spe <- raw_spe


# Create Seurat Object List for PRECAST ----------------------------------


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


# Find gene list from pyschENCODE-spatialDLPFC result

library(tidyverse)
gene_df_raw <- read.csv(
  # TODO: edit the path
  here("code/spatial_clustering/PRECAST",
       "TableS8_sig_genes_FDR5perc_enrichment.csv")
)

# gene_df_raw <- read_csv(
#   here("code/spatial_clustering/PRECAST",
#        "TableS8_sig_genes_FDR5perc_enrichment.csv")
#   )

gene_df <- gene_df_raw |> 
  filter(spatial_domain_resolution == "Sp09") |> 
  group_by(test) |> 
  arrange(fdr, .by_group = TRUE) |> 
  slice_head(n=100)





set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList,
                              selectGenesMethod=NULL,
                              customGenelist=gene_df$ensembl)
# preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, 
                            Sigma_equal = FALSE, coreNum = 8,
                            maxIter = 30, verbose = TRUE)

# K <- as.numeric(Sys.getenv("SGE_TASK_ID"))
K <- 8

# tic()
PRECASTObj <- PRECAST(PRECASTObj, K = seq.int())
# toc()





# Convert PRECAST to a spe object -----------------------------------------

saveRDS(PRECASTObj,
        file = file.path(
          fld_data_spatialcluster, "PRECAST", 
          paste0("test_PRECASTObj_semi_inform_K",K,"_00_model_fitted.rds")
        ))


tmp <- readRDS(file = file.path(
  fld_data_spatialcluster, "PRECAST", 
  paste0("test_PRECASTObj_semi_inform_K",K,"_00_model_fitted.rds")))


# Necessary step to get cluster from resList
PRECASTObj <- SelectModel(PRECASTObj)
# saveRDS(PRECASTObj,
#         file = file.path(
#           fld_data_spatialcluster, "PRECAST", 
#           paste0("test_PRECASTObj_semi_inform_K",K,"_after_SelectModel.rds")
#         ))



# TODO: what does IntegrateSpaData does
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
saveRDS(seuInt,
        file = file.path(
          fld_data_spatialcluster, "PRECAST", 
          paste0("test_seuInt_semi_inform_K",K,"_after_IntegrateSpaData.rds")
        ))

# seuInt <- readRDS(
#   file = file.path(
#     fld_data_spatialcluster, "PRECAST", 
#     paste0("test_seuInt_semi_inform_K",K,"_after_IntegrateSpaData.rds")
#   ))



# Find Marker Genes -------------------------------------------------------

# library(Seurat)
# dat_deg <- FindAllMarkers(seuInt)


# Merge with spe object
stopifnot(ncol(spe) == nrow(seuInt@meta.data))

col_data_df <- seuInt@meta.data |> 
  mutate(cluster = factor(cluster)) |> 
  rename_with(~ paste0("PRECAST_", .x)) |> 
  rownames_to_column(var = "key") |> 
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  ) |> 
  select(-sample_id)

rownames(col_data_df) <- colnames(spe)
# colData(spe) <- NULL
colData(spe) <- DataFrame(col_data_df)



saveRDS(
  spe, 
  file = file.path(
    fld_data_spatialcluster, "PRECAST", 
    paste0("test_PRECASTObj_semi_inform_K",K,".rds")
  )
)


# QC - checking weird clustering samples ----------------------------------
library(ggpubr)
library(escheR)
vis_mito_clus_side_by_side <- function(sampleID, spe){
  ggarrange(
    make_escheR(spe[, spe$sample_id == sampleID]) |> 
      add_fill("PRECAST_cluster"),
    make_escheR(spe[, spe$sample_id == sampleID]) |> 
      add_fill("expr_chrM_ratio"),
    nrow = 1
  )
}

vis_mito_clus_side_by_side(sampleID = "V13F27-293_A1", spe)


# unique(spe$sample_id) |> length()

saveRDS(
  seuInt, 
  file = file.path(
    fld_data_spatialcluster, "PRECAST", 
    paste0("test_seuIntObj_semi_inform_K",K,".rds")
  )
)

# # Visualize Clustering Result ---------------------------------------------
# # TODO: output
library(spatialLIBD)
tmp <- vis_grid_clus(
  spe,
  clustervar = "PRECAST_cluster",
  point_size = 1,
  spatial = FALSE,
  alpha = 1,
  pdf_file = here("plots/spatial_cluster/test_PRECAST_semi_supervised_K8.pdf")
)


vis_grid_clus(
  spe,
  clustervar = "X10x_graphclust",
  point_size = 2,
  spatial = FALSE,
  alpha = 1,
  pdf_file = here("plots/spatial_cluster/triage_10x.pdf")
)


# Session info ------------------------------------------------------------
sessioninfo::session_info()

