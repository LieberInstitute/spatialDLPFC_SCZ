suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(scran)
  library(scuttle)
  library(nnSVG)
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sessioninfo)
})

#TODO: replace this path
spe <- readRDS(here::here(
  "processed-data/rds/spe/spe_after_spot_qc.rds")
)

# Create vector of samples for nnSVG on whole tissue
samples <- unique(spe$sample_id) 


set.seed(12345)
nnsvg_df <- unique(spe$sample_id) |> 
    set_names() |> 
  map_dfc(
    .f = function(sample){
      # browser()
      # TODO: remove 1:5
      spe_sub <- spe[, spe$sample_id == sample]
      
      # TODO: optional
      spe_sub <- filter_genes(spe_sub)
      # re-calculate logcounts after filtering
      spe_sub <- logNormCounts(spe_sub)
      
      # TODO: edit number of thread
      spe_sub <- nnSVG(spe_sub, n_threads = 1)
      
      nnsvg_df <- rowData(spe_sub) |> data.frame() |> 
        select(gene_id, sigma.sq:padj) |> 
        rename_all(~paste0(., "_", sample))
      
      tmp <- rowData(spe) |> data.frame() |> 
        select(gene_id) |> rename_all(~paste0(., "_", sample)) |> 
        left_join(
          nnsvg_df, 
          by = c(paste0("gene_id", "_", sample ))
        )
    }
  )

# Check if merge is correct
stopifnot(all(nnsvg_df$"gene_id_V12F14-053_A1" == nnsvg_df$"gene_id_V12F14-053_B1"))


# save whole tissue nnSVG results
nnsvg_fld <- here::here("processed-data","feature_selection",
                        "nnSVG")

dir.create(nnsvg_fld, showWarnings = FALSE,
           recursive = TRUE)

# Save a data frame that contains only the rowData from nnSVG
# generated spe object
saveRDS(
  nnsvg_df, 
  file = file.path(nnsvg_fld,
                   # TODO: remove test
                  "test_nnSVG.rds")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()