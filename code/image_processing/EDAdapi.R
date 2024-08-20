setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
library("readxl")
library("here")
library("tidyr")
library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("ggridges")
library(escheR)

# spe = readRDS(here("processed-data/rds/spatial_cluster/PRECAST/spe_wo_spg_N63_PRECAST.rds"))
# samples = read_excel(here("raw-data","experiment_info","VisiumSPG_PNN_Master.xlsx"))
# samples <- as.data.frame(samples %>% select('Slide #', 'Array #'))
# samples$sample_id = paste0(samples[,1],"_",samples[,2])
# samples$SPpath = paste0('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/',samples$sample_id,'/outs/spatial/tissue_spot_counts.csv')
# segmentations_list <-
#   lapply(samples$sample_id, function(sample_id) {
#     file <-samples$SPpath[samples$sample_id == sample_id]
#     if (!file.exists(file)) {
#       return(NULL)
#     }
#     x <- read.csv(file)
#     x$key <- paste0(x$barcode, "_", sample_id)
#     return(x)
#   })
#
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <-
#   Reduce(function(...) {
#     merge(..., all = TRUE)
#   }, segmentations_list[lengths(segmentations_list) > 0])
#
#   segmentation_match <- match(spe$key, segmentations$key)
#   segmentation_info <-
#     segmentations[segmentation_match, -which(
#       colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
#     )]
#   colData(spe) <- cbind(colData(spe), segmentation_info)
#
#   saveRDS(spe,here("processed-data", "image_processing", "EDAspe1.rds"))
  
  
  spe = readRDS(here("processed-data", "image_processing", "EDAspe1.rds"))
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
	
  stats_df = read.csv(here("processed-data/image_processing/minmax.csv"))
  coldata_df <- merge(coldata_df, stats_df, by = "sample_id", all.x = TRUE)
  
  coldata_df$iDAPI = coldata_df$IDAPI/coldata_df$DAPI
  coldata_df$iNeuN = coldata_df$INeuN/coldata_df$NeuN
  coldata_df$iWFA = coldata_df$IWFA/coldata_df$WFA
  coldata_df$iClaudin5 = coldata_df$IClaudin5/coldata_df$Claudin5
  
  coldata_df[is.na(coldata_df)] <- 0
  spe$iDAPI = coldata_df$iDAPI
  spe$iNeuN = coldata_df$iNeuN
  spe$iWFA = coldata_df$iWFA
  spe$iClaudin5 = coldata_df$iClaudin5
  
  # Split 'sample_id' at the last underscore and create a new 'slide' column
  coldata_df <- coldata_df %>%
    mutate(slide = sub("_(?!.*_).*", "", sample_id, perl = TRUE))
	
