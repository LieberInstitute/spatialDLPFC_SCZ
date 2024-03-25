# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(RcppML)
  library(Matrix)
  library(scater)
  library(sessioninfo)
  library(tidyverse)
})


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
spe <- readRDS(
  path_PRECAST_int_spe
)

k<- 100

nmf_mdl <- readRDS(
  file = here(
    "processed-data/rds/NMF",
    paste0("test_NMF_SpD5_1-2_k", 100, ".rds")
  )
)

patterns <- t(nmf_mdl$h) # these are the factors
tmp <- rownames_to_column(data.frame(patterns), var = "key")

col_df <- colData(spe) |> data.frame()

fnl_col_df <- left_join(col_df, tmp, by = "key")

rownames(fnl_col_df) |> head()
colData(spe) <- DataFrame(fnl_col_df)
colnames(spe) <- spe$key

.sample_ordered <- metadata(spe)$dx_df |>
  arrange(dx, sample_id) |>
  pull(sample_id)

for (.fac in paste0("NMF_", 1:ncol(patterns))) {
  # .fac <- "NMF_1"
  vis_grid_gene(
    spe,
    geneid = .fac,
    pdf_file = here(paste0("plots/NMF/SpD5/test_spot_plot_SpD_", .fac, ".pdf")),
    sample_order = .sample_ordered,
    point_size = 0.8,
    spatial = FALSE
  )
}

