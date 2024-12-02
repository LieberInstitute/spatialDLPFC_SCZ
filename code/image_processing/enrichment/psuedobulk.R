setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/')
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_w_spg_N63.rds"))

## Load SpD data ----
finalized_spd <- readRDS(here("processed-data/rds/spatial_cluster", "PRECAST","test_clus_label_df_semi_inform_k_2-16.rds"))

## Attach SpD label to spe ----
col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(spe$spg_PDAPI > 0.05,FALSE, TRUE)
spe$neun_pos <- ifelse(spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,TRUE, FALSE)
spe$vasc_pos <- ifelse(spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,TRUE, FALSE)

spe_ntc <- spe[, colData(spe)$dx == "ntc"]
head(colData(spe_ntc))


 ## Per sample & domain ----
neuropil_pseudo <-
  registration_pseudobulk(
    spe_ntc,
    var_registration = "neuropil_pos",
    var_sample_id = "sample_id",
    covars = c("age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here("processed-data", "image_processing", "enrichment", "neuropil_pseudo.rds")
    )

dim(neuropil_pseudo)
#[1] 15717    62	

neun_pseudo <-
  registration_pseudobulk(
    spe_ntc,
    var_registration = "neun_pos",
    var_sample_id = "sample_id",
    covars = c("age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here("processed-data", "image_processing", "enrichment", "neun_pseudo.rds")
    )

dim(neun_pseudo)
#[1] 15637    62	

pnn_pseudo <-
  registration_pseudobulk(
    spe_ntc,
    var_registration = "pnn_pos",
    var_sample_id = "sample_id",
    covars = c("age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here("processed-data", "image_processing", "enrichment", "pnn_pseudo.rds")
    )

dim(pnn_pseudo)
#[1] 14899    62	

vasc_pseudo <-
  registration_pseudobulk(
    spe_ntc,
    var_registration = "vasc_pos",
    var_sample_id = "sample_id",
    covars = c("age", "sex", "lot_num", "slide_id"),
    min_ncells = 10,
    pseudobulk_rds_file = here("processed-data", "image_processing", "enrichment", "vasc_pseudo.rds")
    )

dim(vasc_pseudo)
#[1] 14755    62	

