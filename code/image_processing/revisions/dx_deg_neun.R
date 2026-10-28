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
raw_spe <- readRDS(here("processed-data/rds/01_build_spe", "fnl_spe_kept_spots_only.rds"))
colData(raw_spe) |> names()

## Subset to neuropil+ only object ----
spg_spe <- raw_spe[, raw_spe$neun_pos &!raw_spe$vasc_pos &!raw_spe$neuropil_pos]

# gene by spots
# dim: 36601 55017 

#dim(spg_spe)
#[1] 36601 38129
# Create pseudobulked data ----
sce_pseudo <-
  registration_pseudobulk(
    spg_spe,
    var_registration = "fnl_spd",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10)

pb_neun <- sce_pseudo[,sce_pseudo$registration_variable %in%sprintf("spd%02d", c(2, 3, 5, 6, 7))]

ncol(pb_neun)
# [1] 251

# DE analysis ----
dx_mod <-
  registration_model(
    pb_neun,
    covars = c("age", "sex"),
    var_registration = "dx"
  )

# 2025-04-14 11:45:46.326476 create model matrix

dx_block_cor <- registration_block_cor(
  pb_neun,
  registration_model = dx_mod,
  var_sample_id = "sample_id"
)

# 2025-04-14 11:45:51.214688 run duplicateCorrelation()

dx_res <- registration_stats_enrichment(
  pb_neun,
  block_cor = dx_block_cor,
  covars = c("age", "sex", "PRECAST_07"),
  var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

dx_res |> arrange(fdr_scz) |> head(10)
dx_res |> arrange(t_stat_scz) |> head(10)
dx_res |> arrange(desc(t_stat_scz)) |> head(10)
dx_res |> arrange(desc(logFC_scz)) |> head(10)

# 2025-04-14 11:47:44.999823 computing enrichment statistics
# 2025-04-14 11:47:46.405967 extract and reformat enrichment results

## Save results ----
write_csv(dx_res |> arrange(fdr_scz),here("processed-data/image_processing/revisions","neun-dx_DEG-GM.csv"))
