# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  # library(limma)
  library(sessioninfo)
  library(here)
})

# Load Data ----
## Load SPG spe object ----
raw_spe <- readRDS(
  here::here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

colData(raw_spe) |> names()

## Subset to vasc+ only object ----
spg_spe <- raw_spe[, raw_spe$vasc_pos == TRUE]
# gene by spots
# dim: 36601 21906

# Create pseudobulked data ----
sce_pseudo <-
  registration_pseudobulk(
    spg_spe,
    var_registration = "fnl_spd",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    # NOTE: not named pseudo_vasc_pos_donor_spd.rds, the 2024 object still read by
    # 15_eqtl_coloc and code_for_isee_app. Same counts and samples as this object,
    # but lacks some colData (e.g. slide_id, fnl_spd)
    pseudobulk_rds_file = here(
      "processed-data/rds/PB_dx_spg",
      "pseudo_vasc_pos_donor_spd_fnl.rds"
    )
  )



# Session Info ----
session_info()