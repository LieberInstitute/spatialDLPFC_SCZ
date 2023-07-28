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
  "processed-data/rds/spe/01_build_spe/spe_raw.rds")
)

# Convert to seuList
seuList <- (unique(spe$sample_id) |> 
  # TODO: remove the subsetting
  set_names(unique(spe$sample_id)))[1:2] |> 
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
gene_df_raw <- read_csv(here("code/visium_spatial_clustering/PRECAST/TableS8_sig_genes_FDR5perc_enrichment.csv"))

gene_df <- gene_df |> 
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

K <- as.numeric(Sys.getenv("SGE_TASK_ID"))

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

save(PRECASTObj, file = here("processed-data", "06_clustering", "PRECAST", paste0("allSamples_PRECASTObj_nnSVG_2500_",K,".Rdata")))
