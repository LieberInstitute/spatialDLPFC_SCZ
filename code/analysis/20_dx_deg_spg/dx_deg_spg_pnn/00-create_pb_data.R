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

## Subset to neuropil+ only object ----
spg_spe <- raw_spe[, raw_spe$pnn_pos == TRUE]
# gene by spots
# dim: 36601 55017 

# Create pseudobulked data ----
sce_pseudo <-
  registration_pseudobulk(
    spg_spe,
    var_registration = "fnl_spd",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here(
      "processed-data/rds/PB_dx_spg",
      "pseudo_pnn_pos_donor_spd.rds"
    )
  )

# 2025-07-02 16:03:08.208829 make pseudobulk object
# 2025-07-02 16:03:10.438186 dropping 147 pseudo-bulked samples that are below 'min_ncells'.
# 2025-07-02 16:03:10.467088 drop lowly expressed genes
# 2025-07-02 16:03:10.588776 normalize expression
# 2025-07-02 16:03:11.478425 saving sce_pseudo to /Users/bguo6/GitHub/spatialDLPFC_SCZ/processed-data/rds/PB_dx_spg/pseudo_pnn_pos_donor_spd.rds


# Session Info ----
session_info()