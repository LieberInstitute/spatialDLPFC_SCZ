setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  # library(limma)
  library(sessioninfo)
  library(here)
})

## Load SPG spe object ----
raw_spe <- readRDS(here("processed-data/rds/01_build_spe","fnl_spe_kept_spots_only.rds"))
colData(raw_spe) |> names()

## Subset to pnn+ only object ----
spg_spe <- raw_spe[, raw_spe$pnn_pos == TRUE]
# dim: 36601 55017 , [1] 36601 12333

spg_spe <- raw_spe[, raw_spe$pnn_pos &!raw_spe$neuropil_pos &!raw_spe$vasc_pos]
# dim: 36601  6347

# Create pseudobulked data ----
sce_pseudo <-
  registration_pseudobulk(
    spg_spe,
    var_registration = "fnl_spd",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10
  )

#  2026-09-29 13:21:16.785947 make pseudobulk object
#  2026-09-29 13:21:19.546606 dropping 150 pseudo-bulked samples that are below 'min_ncells'.
#  2026-09-29 13:21:19.617636 drop lowly expressed genes
#  2026-09-29 13:21:19.83495 normalize expression

pb_pnn <- sce_pseudo[,sce_pseudo$registration_variable %in%sprintf("spd%02d", c(2, 3, 5, 6, 7))]
ncol(pb_pnn)
#[1] 159

dx_mod <-registration_model(pb_pnn,covars = c("age", "sex"),var_registration = "dx")
dx_block_cor <- registration_block_cor(pb_pnn,registration_model = dx_mod,var_sample_id = "sample_id")
dx_res <- registration_stats_enrichment(
  pb_pnn,
  block_cor = dx_block_cor,
  covars = c("age", "sex", "PRECAST_07"),
  var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

write_csv(dx_res |> arrange(fdr_scz),here("processed-data/image_processing/revisions","pnn-dx_DEG-GM.csv"))