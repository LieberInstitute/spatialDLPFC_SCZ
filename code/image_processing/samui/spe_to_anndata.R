setwd("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/")
suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

spe <- readRDS(here("processed-data", "image_processing", "EDAspe.rds"))
finalized_spd <- readRDS(here("processed-data/rds/spatial_cluster","PRECAST","test_clus_label_df_semi_inform_k_2-16.rds"))
coldata_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )
 rownames(coldata_df) <- colnames(spe)
 colData(spe) <- DataFrame(coldata_df)
 
 df = as.data.frame(colData(spe))
 df <- df %>% mutate(PNN = ifelse(PWFA > 0.05, "PNN+", "PNN-"))
 df <- df %>% mutate(DAPI = ifelse(PDAPI > 0.05 & PDAPI < 0.5 , "DAPI+", "DAPI-"))
 df <- df %>% mutate(NeuN = ifelse(PNeuN > 0.05 & PNeuN < 0.3 , "NeuN+", "NeuN-"))
 df <- df %>% mutate(Claudin = ifelse(PClaudin5 > 0.05 & PClaudin5 < 0.20, "Claudin+", "Claudin-"))
 
 colData(spe) <- DataFrame(df)
 
spe_out <- here("processed-data", "image_processing", "spg.h5ad")

write_anndata <- function(sce, out_path) {
  invisible(
    basiliskRun(
      fun = function(sce, filename) {
        library("zellkonverter")
        library("reticulate")
        
        # Convert SCE to AnnData:
        adata <- SCE2AnnData(sce)
        
        #  Write AnnData object to disk
        adata$write(filename = filename)
        
        return()
      },
      env = zellkonverterAnnDataEnv(),
      sce = sce,
      filename = out_path
    )
  )
}


#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnData, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe)$spatial <- spatialCoords(spe)

write_anndata(spe, spe_out)

session_info()
