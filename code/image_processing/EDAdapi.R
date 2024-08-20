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
	
    ##DAPI
    p1 = ggplot(coldata_df, aes(x = PDAPI)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
    p11 = ggplot(coldata_df, aes(x = iDAPI)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
	  
    p2 = ggplot(coldata_df, aes(x = PDAPI)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
    p22 = ggplot(coldata_df, aes(x = iDAPI)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
		  
    p3 = ggplot(coldata_df, aes(x = PDAPI,y=iDAPI, color = sample_id)) + theme_minimal() + geom_point(alpha = 0.3) +theme(legend.position = "none")
    p33 = ggplot(coldata_df, aes(x = PDAPI,y=iDAPI)) + geom_point(color = "black") + theme_minimal() + scale_y_log10() + scale_x_log10()
		  
    ggsave(here("plots", "image_processing", "DAPIp1.png"), plot = p1, width = 6, height = 4, dpi = 300)
    ggsave(here("plots", "image_processing", "DAPIp11.png"), plot = p11, width = 6, height = 4, dpi = 300)
    ggsave(here("plots", "image_processing", "DAPIp2.png"), plot = p2, width = 6, height = 4, dpi = 300)
    ggsave(here("plots", "image_processing", "DAPIp22.png"), plot = p22, width = 6, height = 4, dpi = 300)
    ggsave(here("plots", "image_processing", "DAPIp3.png"), plot = p3, width = 6, height = 4, dpi = 300)
    ggsave(here("plots", "image_processing", "DAPIp33.png"), plot = p33, width = 6, height = 4, dpi = 300)
  
 