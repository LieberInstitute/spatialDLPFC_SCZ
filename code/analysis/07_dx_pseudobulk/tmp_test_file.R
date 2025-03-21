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
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

spe$in_tissue |> table()

paste0(spe$sample_id, ".", spe$PRECAST_07) |>
  unique()




sum(is.na(spe$PRECAST_07))

df_miss_PRECAST <- colData(spe) |> data.frame() |>

df_miss_PRECAST$sample_id |> unique()
df_miss_PRECAST$sample_label |> unique()

df_miss_PRECAST$sample_label |> table()

colData(spe) |> data.frame() |>
  filter(sample_id == "V13M06-344_B1") |> 
  select(PRECAST_07) |> View()



## Create donor-spd pseudobulk data ----
sce_pseudo_tmp <-
  registration_pseudobulk(
    spe,
    var_registration = "spd_label",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 0#,
    # pseudobulk_rds_file = here(
    #   "processed-data/rds/07_dx_pseudobulk",
    #   "sce_pseudo_PRECAST07_donor_spd.rds"
    # )
  )

tmp_df <- colData(sce_pseudo_tmp) |> data.frame() |>
transmute(
  pb_id  = paste0(sample_id, ".", PRECAST_07),
  ncells
)

 colData(sce_pseudo_tmp) |> data.frame() |> group_by(PRECAST_07) |> summarise(n = n())
