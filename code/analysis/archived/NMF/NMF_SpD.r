# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(RcppML)
  library(Matrix)
  library(scater)
  library(sessioninfo)
  library(tidyverse)
})


# Config -----------------------------------------------------------
## Path Config -------------------------------------------------------------
# TODO: decide if this needs to be revised
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster"
)

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform", ".rds")
)

# Load data ---------------------------------------------------------------
raw_spe <- readRDS(
  path_PRECAST_int_spe
)

spe <- raw_spe[, raw_spe$PRECAST_5 %in% c(1, 2)]
# spe <- raw_spe[, raw_spe$PRECAST_5 %in% c(1, 2) & raw_spe$sample_id == "V12F14-053_A1"] # Test ony



# NOTE: For testing only
# spe_backup <- spe
# subset_id <- metadata(spe)$dx_df |> group_by(dx) |>
#   slice_head(n=2) |> ungroup() |>
#   pull(sample_id)
#
# spe <- spe_backup[ ,spe_backup$sample_id %in% subset_id]

## Add dx informaiton to the dataset
# spe$dx <- metadata(spe)$dx_df$dx[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )]
#
# spe$age <- metadata(spe)$dx_df$age[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )]
#
# spe$sex <- metadata(spe)$dx_df$sex[
#   match(
#     spe$sample_id,
#     metadata(spe)$dx_df$sample_id
#   )]



# Subsetting one dx group -------------------------------------------------
# spe <- spe[, spe$dx == unique(spe$dx)[1]]

cat("Total number of spots is ", ncol(spe), ".\n")
# NMF ---------------------------------------------------------------------

# library(spatialLIBD)
# library(escheR)
# library(here)
# library(tidyverse)
# library(ggforce)

# Code adapted from Cindy's code
# https://github.com/cindyfang70/xenium-sandbox/blob/main/code/NMF/01_exploratory_visium_nmf.R


## Config ------------------------------------------------------------------
k <- 100
spe <- spe[, spe$sum_umi !=0]
spe <- logNormCounts(spe)

A <- logcounts(spe) # using logcounts because there are multiple datasets

model <- RcppML::nmf(A, k = k, seed = 1237)

# Assign key
# model$h |> dim()
rownames(model$h) <- paste0("NMF_", seq.int(nrow(model$h)))
colnames(model$h) <- spe$key


## Save results
saveRDS(
  model,
  file = here(
    "processed-data/rds/NMF",
    paste0("test_NMF_SpD5_1-2_k", k, ".rds")
  )
)

print("Finish NMF")


# Session Info ------------------------------------------------------------

sessioninfo::session_info()
