# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(sessioninfo)
})


# Path --------------------------------------------------------------------
#TODO: doesn't seem to use this. Delete
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

# Load QC-ed SPE object --------------------------------------------------
raw_spe <- readRDS(here::here(
  "processed-data/rds/",
  #TODO: replace this path
  "test_spe_after_spot_qc_63.rds")
)

spe <- raw_spe[, raw_spe$sample_id %in% c("V13F27-294_B1", "V12F14-057_A1")]

library(scuttle)
# Normalization
spe <- logNormCounts(spe, log = FALSE)

stopifnot("normcounts" %in% assayNames(spe))
# Investigation -----------------------------------------------------------


library(spatialLIBD)
### Mito_Ratio --------------------------------------------------------------
vis_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",    # PCP4
  # geneid = "spg_NDAPI",
  spatial=FALSE,
  assayname = "normcounts"
)

### WM (MOBP) --------------------------------------------------------------
# rowData(spe)["ENSG00000168314",]
vis_gene(
  spe = spe,
  geneid = "ENSG00000168314",    # MOBP
  # geneid = "spg_NDAPI",
  spatial=FALSE,
  assayname = "normcounts"
)


### Gray Matter (SNAP25) -----------------------------------------------------------
vis_gene(
  spe = spe, sampleid = "V12F14-057_A1", 
  geneid = "ENSG00000132639",    # SNAP25
  # geneid = "spg_NDAPI",
  spatial=FALSE,
  assayname = "normcounts"
)



### (Layer 1) AQP4 --------------------------------------------------------------
library(ggpubr)

ggarrange(
  vis_gene(
    spe = spe, sampleid = "V13F27-294_B1",
    geneid = "ENSG00000171885",    # AQP4
    # geneid = "spg_NDAPI",
    spatial=FALSE,
    assayname = "normcounts"
  ),
  vis_gene(
    spe = spe, sampleid = "V12F14-057_A1",
    geneid = "ENSG00000171885",    # AQP4
    # geneid = "spg_NDAPI",
    spatial=FALSE,
    assayname = "normcounts"
  ), 
  nrow = 1
)


### PCP4 --------------------------------------------------------------

ggarrange(
  vis_gene(
    spe = spe, sampleid = "V13F27-294_B1",
    geneid = "ENSG00000183036",     # PCP4
    # geneid = "spg_NDAPI",
    spatial=FALSE,
    assayname = "normcounts"
  ),
  vis_gene(
    spe = spe, sampleid = "V12F14-057_A1",
    geneid = "ENSG00000183036",     # PCP4
    # geneid = "spg_NDAPI",
    spatial=FALSE,
    assayname = "normcounts"
  ), 
  nrow = 1
)




# Session info ------------------------------------------------------------
sessioninfo::session_info()
