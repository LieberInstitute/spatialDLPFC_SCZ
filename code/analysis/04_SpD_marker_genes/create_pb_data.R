# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(sessioninfo)
})

# Load data -----
## Load spe object ----
spe <- readRDS(
  here::here(
    "processed-data/rds/02_visium_qc",
    "test_qc_spe_w_spg_N63_no_img.rds"
  )
)

PRECAST_df <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Merge Demo info ----
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$age <- metadata(spe)$dx_df$age[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$slide_id <- sapply(strsplit(spe$sample_id, "_"), function(x) x[1])
spe$lot_num <- sapply(strsplit(spe$slide_id, "-"), function(x) x[1])

## Merge SPD clus ----
stopifnot(nrow(PRECAST_df) == ncol(spe))

precast_vars <- grep(
  "^PRECAST_", colnames(PRECAST_df),
  value = TRUE
)
spe <- spe[, spe$key %in% PRECAST_df$key]
# raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
col_data_df <- PRECAST_df |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )
rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)


# Enrichment test for each SpD ----
vars <- grep("^PRECAST", names(colData(spe)), value = TRUE)

# .var <- vars[7]

for (.var in vars) {
  cat(
    "Start building pseudobulk makers for ",
    .var, "\n"
  )
  stopifnot(
    is.character(spe[[.var]])
  )
  
  sce_pseudo <-
    registration_pseudobulk(
      spe,
      var_registration = .var,
      var_sample_id = "sample_id",
      covars = c("dx", "age", "sex", "lot_num", "slide_id"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "rds", "layer_spd",
        paste0("test_spe_pseudo_", .var, ".rds")
      )
    )
}

# Session Info ----
sessioninfo::session_info()
