setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  library(limma)
  library(sessioninfo)
  library(here)
  library(ggrepel)
})

## Load SPG spe object ----
raw_spe <- readRDS(here("processed-data/rds/01_build_spe","fnl_spe_kept_spots_only.rds"))
colData(raw_spe) |> names()

## Subset to neuropil+ only object ----
spg_spe <- raw_spe[, raw_spe$vasc_pos == TRUE]
# dim: [1] 36601 21906 

spg_spe <- raw_spe[, raw_spe$vasc_pos &!raw_spe$neun_pos &!raw_spe$pnn_pos &!raw_spe$neuropil_pos]
# dim: 36601 12750

# Create pseudobulked data ----
sce_pseudo <-
  registration_pseudobulk(
    spg_spe,
    var_registration = "fnl_spd",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10)
	
#2026-09-29 12:42:09.431836 make pseudobulk object
#2026-09-29 12:42:18.720315 dropping 8 pseudo-bulked samples that are below 'min_ncells'.
#2026-09-29 12:42:18.903045 drop lowly expressed genes
#2026-09-29 12:42:19.484341 normalize expression

pb_vasc <- sce_pseudo[,sce_pseudo$registration_variable %in%sprintf("spd%02d", c(2, 3, 5, 6, 7))]
ncol(pb_vasc)
ncol(sce_pseudo)
#[1] 314

dx_mod <-registration_model(sce_pseudo,covars = c("age", "sex"), var_registration = "dx")
dx_block_cor <- registration_block_cor(sce_pseudo, registration_model = dx_mod, var_sample_id = "sample_id")
dx_res <- registration_stats_enrichment(sce_pseudo,
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

write_csv(dx_res |> arrange(fdr_scz),here("processed-data/image_processing/revisions","vasc-dx_DEG.csv"))