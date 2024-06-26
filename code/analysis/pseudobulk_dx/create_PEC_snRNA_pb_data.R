library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

data_dir <- file.path(
  "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/",
  "processed-data", "rdata", "spe", "12_spatial_registration_sn"
)

#### Load sn data & exclude drop cells ####
load(file = here(data_dir, "sce_DLPFC.Rdata"), verbose = TRUE)
sce <- sce[, sce$cellType_hc != "drop"]
sce$cellType_hc <- droplevels(sce$cellType_hc)

table(sce$cellType_hc)

## Factor categorical variables used as covariates
colData(sce)$Position <- as.factor(colData(sce)$Position)
colData(sce)$sex <- as.factor(colData(sce)$sex)

## Use all unique ensembl IDs as rownames
rownames(sce) <- rowData(sce)$gene_id

sce_pb <- registration_pseudobulk(
  sce = sce,
  var_registration = "cellType_hc",
  var_sample_id = "Sample",
  covars = c("Position", "age", "sex"),
  # pseudobulk_rds_file = here(
  #   "code/analysis/pseudobulk_dx",
  #   "spatialDLPFC_sn_hc_pb_data.rds"
  # )
)
