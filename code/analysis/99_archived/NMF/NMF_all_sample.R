# Load Libray -----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(scuttle)
  library(RcppML)
  library(Matrix)
  library(sessioninfo)
})

# Load data -----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)

# NMF -----
# Code adapted from Cindy Fang's code
# https://github.com/cindyfang70/xenium-sandbox/blob/main/code/NMF/01_exploratory_visium_nmf.R


## Config ----
stopifnot("logcounts" %in% assayNames(spe))

k <- 100
A <- logcounts(spe)

stopifnot(sum(A<0) <= 0)

model <- RcppML::nmf(A, k = k, seed = 20240516)
print("Finish NMF")

## Save results
saveRDS(
  model,
  file = here(
    "processed-data/rds/NMF",
    paste0("test_NMF_all_k", k, ".rds")
  )
)

print("Finish save NMF")


# Session Info -----
sessioninfo::session_info()
