# Load packages----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
})

# Load Data ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)

# TODO: read in SPD data ----


# Create channel specific PB ----

# TODO: add an for loop to iterate over each SPG chanel
spg_names <- c("pnn_pos", "neuropil_pos", "neun_pos", "vasc_pos")

## Create positive PNN only spe ----
for (.spg in spg_names) {
  .spg <- spg_names[1]

  pnn_spe <- spe[, spe[[.spg]] == TRUE]

  ## Pseudobulk based on individual
  sce_pseudo <-
    registration_pseudobulk(
      pnn_spe,
      var_registration = "pnn_pos",
      var_sample_id = "sample_id",
      covars = c("dx", "age", "sex", "lot_num", "slide_id"),
      min_ncells = 10
      # TODO: change path
      # pseudobulk_rds_file = here(
      #   "processed-data", "rds", "layer_spd",
      #   paste0("test_spe_pseudo_", .var, ".rds")
      # )
    )

  ## Per sample & domain ----
  # TODO: recall how to create pseudobulked data by sample_id and SpD
  # TODO: save pseudobulked data
  sce_pseudo <-
    registration_pseudobulk(
      pnn_spe,
      var_registration = "pnn_pos",
      # TODO: add spd to pseudobulk
      var_sample_id = "sample_id",
      covars = c("dx", "age", "sex", "lot_num", "slide_id"),
      min_ncells = 10
      # TODO: change path
      # pseudobulk_rds_file = here(
      #   "processed-data", "rds", "layer_spd",
      #   paste0("test_spe_pseudo_", .var, ".rds")
      # )
    )
}


# Session Info ----
sessioninfo::session_info()
