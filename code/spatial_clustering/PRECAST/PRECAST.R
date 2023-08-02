library(here)
library("dplyr")
library("purrr")
library("Seurat")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("PRECAST")
# library("tictoc")



#TODO: replace this path
spe <- readRDS(here::here(
  "processed-data/rds/spe/spe_after_spot_qc.rds")
)

# Convert to seuList
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


# Find gene list from pyschENCODE-spaitalDLPFC result

library(tidyverse)
gene_df_raw <- read_csv(here("code/spatial_clustering/PRECAST/TableS8_sig_genes_FDR5perc_enrichment.csv"))

gene_df <- gene_df_raw |> 
  filter(spatial_domain_resolution == "Sp09") |> 
  group_by(test) |> 
  arrange(fdr, .by_group = TRUE) |> 
  slice_head(n=100)




  
set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList,
                              selectGenesMethod=NULL,
                              customGenelist=gene_df$ensembl)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 8, maxIter = 30, verbose = TRUE)

# K <- as.numeric(Sys.getenv("SGE_TASK_ID"))
K <- 8

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

dir.create(
  here("processed-data", "clustering", "PRECAST"),
  recursive = T, showWarnings = FALSE
)




# Convert PRECAST to a spe object -----------------------------------------

# Necessary step to get cluster from resList
PRECASTObj <- SelectModel(PRECASTObj)

seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
col_data_df <- seuInt@meta.data |> 
  mutate(cluster = factor(cluster)) |> 
  rename_with(~ paste0("PRECAST_", .x)) |> 
  rownames_to_column(var = "key") |> 
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)

colData(spe) <- DataFrame(col_data_df)

saveRDS(spe, file = here("processed-data", "clustering", "PRECAST", 
                         paste0("test_PRECASTObj_semi_inform_K",K,".rds")))


# Visualize Clustering Result ---------------------------------------------
# TODO: output
tmp <- vis_grid_clus(spe,
              clustervar = "PRECAST_cluster", return_plots = TRUE)


