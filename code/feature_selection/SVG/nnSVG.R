suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(scran)
  library(nnSVG)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sessioninfo)
})

#TODO: replace this path
spe <- readRDS(here::here(
  "processed-data/rds/spe/01_build_spe/spe_raw.rds")
)

# Create vector of samples for nnSVG on whole tissue
samples <- unique(spe$sample_id)
# Run nnSVG once per sample whole tissue and store lists of top SVGs
res_list <- as.list(rep(NA, length(samples)))

names(res_list) <- samples



for (s in samples) {
  
  # select sample_id
  ix <- spe$sample_id == s
  spe_sub <- spe[, ix]
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(spe_sub)
  
  # re-calculate logcounts after filtering
  spe_sub <- logNormCounts(spe_sub)
  
  # run nnSVG
  set.seed(12345)
  message('running nnSVG')
  spe_sub <- nnSVG(spe_sub, n_threads = 10)
  
  # store whole tissue results
  message('saving data')
  res_list[[s]] <- rowData(spe_sub)
  # temp = rowData(spe_sub)
  # save(temp, file = here::here("processed-data","07_Feature_selection", "nnSVG", paste0(samples[s], ".Rdata")))
}

dir.create(
  here("processed-data","feature_selection", "nnSVG"),
  recursive = T, showWarnings = FALSE
)


# save whole tissue nnSVG results
saveRDS(res_list, file = here::here("processed-data","feature_selection", "nnSVG", "nnSVG.rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()