setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
library("readxl")
library("here")
library("tidyr")
library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("ggridges")
library(escheR)
library("purrr")
library(tidyverse)
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
#   saveRDS(spe,here("processed-data", "image_processing", "EDAspe_clean.rds"))
# 
  
  spe = readRDS(here("processed-data", "image_processing", "EDAspe_clean.rds"))
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
  
    filtered_data <- coldata_df[coldata_df$iDAPI < 0.12 | coldata_df$PDAPI > 0.75 , ]	

    p3 = ggplot(coldata_df, aes(x = PDAPI,y=iDAPI, color = slide)) + geom_point(alpha = 0.3) +labs(shape = "Slide") +
         geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +  # Horizontal line at y = 0.01
         geom_vline(xintercept = 0.03, linetype = "dashed", color = "red") +
  	   geom_text(data = filtered_data, aes(label = sample_id), vjust = -1, color = "black") + theme_minimal() # Vertical line at x = 0.05

   ggsave(here("plots", "image_processing", "DAPIsampleexplore.png"), plot = p3, width = 12, height = 4, dpi = 300)
 
   library(escheR)
   spd_anno_df <- read_csv(here("processed-data/man_anno","spd_labels_k7.csv")) |>
     mutate(anno_lab = paste0(label, " (", spd, ") "))
  
   spe$dapi_pos <- spe$PDAPI >= 0.75  
   spe$dapi_pos <- spe$iDAPI <= 0.25
   sample_id = "V13M06-280_A1"
   plot(
      make_escheR(spe[, spe$sample_id == sample_id]) |>
        add_ground("PRECAST_07") |>
        add_fill("dapi_pos") +
        scale_color_manual(
          values = set_names(
            Polychrome::palette36.colors(7)[seq.int(7)],
            unique(spe[["PRECAST_07"]]) |> sort()
          ),
          limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
          labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
        ) +
        scale_fill_manual(
          values = c(
            "TRUE" = "black",
            "FALSE" = "transparent"
          )
        ) +
        labs(title = sample_id) )
	  
	    ##NeuN
	    p1 = ggplot(coldata_df, aes(x = PNeuN)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
	    p11 = ggplot(coldata_df, aes(x = iNeuN)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
	  
	    p2 = ggplot(coldata_df, aes(x = PNeuN)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
	    p22 = ggplot(coldata_df, aes(x = iNeuN)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
		  
	    p3 = ggplot(coldata_df, aes(x = PNeuN,y=iNeuN, color = sample_id)) + geom_point(alpha = 0.3) + theme_minimal() +theme(legend.position = "none")  
	    p33 = ggplot(coldata_df, aes(x = PNeuN,y=iNeuN)) + geom_point(color = "black") + theme_minimal() + scale_y_log10() + scale_x_log10()
		  
	    ggsave(here("plots", "image_processing", "NeuNp1.png"), plot = p1, width = 6, height = 4, dpi = 300)
	    ggsave(here("plots", "image_processing", "NeuNp11.png"), plot = p11, width = 6, height = 4, dpi = 300)
	    ggsave(here("plots", "image_processing", "NeuNp2.png"), plot = p2, width = 6, height = 4, dpi = 300)
	    ggsave(here("plots", "image_processing", "NeuNp22.png"), plot = p22, width = 6, height = 4, dpi = 300)
	    ggsave(here("plots", "image_processing", "NeuNp3.png"), plot = p3, width = 6, height = 4, dpi = 300)
	    ggsave(here("plots", "image_processing", "NeuNp33.png"), plot = p33, width = 6, height = 4, dpi = 300)
	
	    filtered_data <- coldata_df[coldata_df$PNeuN > 0.75 & coldata_df$iNeuN < 0.15, ]	
	     p3 = ggplot(coldata_df, aes(x = PNeuN,y=iNeuN, color = slide)) + geom_point(alpha = 0.3) +labs(shape = "Slide") +
	          geom_hline(yintercept = 0.12, linetype = "dashed", color = "red") +  # Horizontal line at y = 0.01
	          geom_vline(xintercept = 0.03, linetype = "dashed", color = "red") +
	   	   geom_text(data = filtered_data, aes(label = sample_id), vjust = -1, color = "black") + theme_minimal() # Vertical line at x = 0.05

	    ggsave(here("plots", "image_processing", "NeuNsampleexplore.png"), plot = p3, width = 12, height = 4, dpi = 300)
   
	    spe$dapi_pos <- spe$PClaudin5 >= 0.1  
	    #spe$dapi_pos <- spe$iDAPI <= 0.25
	    sample_id = "V12D07-334_B1"
	    plot(
	       make_escheR(spe[, spe$sample_id == sample_id]) |>
	         add_ground("PRECAST_07") |>
	         add_fill("dapi_pos") +
	         scale_color_manual(
	           values = set_names(
	             Polychrome::palette36.colors(7)[seq.int(7)],
	             unique(spe[["PRECAST_07"]]) |> sort()
	           ),
	           limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
	           labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
	         ) +
	         scale_fill_manual(
	           values = c(
	             "TRUE" = "black",
	             "FALSE" = "transparent"
	           )
	         ) +
	         labs(title = sample_id) )
	  
   
 
		     ##WFA
		     p1 = ggplot(coldata_df, aes(x = PWFA)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
		     p11 = ggplot(coldata_df, aes(x = iWFA)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
	  
		     p2 = ggplot(coldata_df, aes(x = PWFA)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
		     p22 = ggplot(coldata_df, aes(x = iWFA)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
   
		     p3 = ggplot(coldata_df, aes(x = PWFA,y=iWFA, color = sample_id)) + geom_point(alpha = 0.3)+ theme_minimal() +theme(legend.position = "none") 
		     p33 = ggplot(coldata_df, aes(x = PWFA,y=iWFA)) + geom_point(color = "black") + theme_minimal() + scale_y_log10() + scale_x_log10()
		  
		     ggsave(here("plots", "image_processing", "WFAp1.png"), plot = p1, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "WFAp11.png"), plot = p11, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "WFAp2.png"), plot = p2, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "WFAp22.png"), plot = p22, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "WFAp3.png"), plot = p3, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "WFAp33.png"), plot = p33, width = 6, height = 4, dpi = 300)

		     ##CLaudin
		     p1 = ggplot(coldata_df, aes(x = PClaudin5)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
		     p11 = ggplot(coldata_df, aes(x = iClaudin5)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()
	  
		     p2 = ggplot(coldata_df, aes(x = PClaudin5)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
		     p22 = ggplot(coldata_df, aes(x = iClaudin5)) + geom_histogram(binwidth = 0.001, color = "black") +theme_minimal()+ scale_y_log10() 
		  
		     p3 = ggplot(coldata_df, aes(x = PClaudin5,y=iClaudin5, color = sample_id)) + geom_point(alpha = 0.3) + theme_minimal() +theme(legend.position = "none") 
		     p33 = ggplot(coldata_df, aes(x = PClaudin5,y=iClaudin5)) + geom_point(color = "black") + theme_minimal() + scale_y_log10() + scale_x_log10()
		  
		     ggsave(here("plots", "image_processing", "Claudin5p1.png"), plot = p1, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "Claudin5p11.png"), plot = p11, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "Claudin5p2.png"), plot = p2, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "Claudin5p22.png"), plot = p22, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "Claudin5p3.png"), plot = p3, width = 6, height = 4, dpi = 300)
		     ggsave(here("plots", "image_processing", "Claudin5p33.png"), plot = p33, width = 6, height = 4, dpi = 300)

 
 

		     filtered_data <- coldata_df[coldata_df$PWFA > 0.25 & coldata_df$iWFA <=0.37, ]	
			 filtered_data <- coldata_df[coldata_df$P > 0.25 & coldata_df$iWFA <=0.37, ]	
		     #filtered_data <- coldata_df[coldata_df$sample_id == "V13F27-296_B1", ]

		     p3 = ggplot(coldata_df, aes(x = PClaudin5,y=iClaudin5, color = sample_id)) + geom_point(alpha = 0.3) +#labs(shape = "Slide") +
		          geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +  # Horizontal line at y = 0.01
		          geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
		   	   geom_text(data = coldata_df, aes(label = sample_id), vjust = -1, color = "black") + theme_minimal()+
			   xlim(c(0.1,0.5)) + ylim(c(0,0.1))+theme(legend.position = "none")  # Vertical line at x = 0.05

		   	   spe$dapi_pos <- spe$PWFA >= 0.05 
		   	   #spe$dapi_pos <- spe$iDAPI <= 0.25
		   	   sample_id = "V12F14-053_A1"
		   	   plot(
		   	      make_escheR(spe[, spe$sample_id == sample_id]) |>
		   	        add_ground("PRECAST_07") |>
		   	        add_fill("dapi_pos") +
		   	        scale_color_manual(
		   	          values = set_names(
		   	            Polychrome::palette36.colors(7)[seq.int(7)],
		   	            unique(spe[["PRECAST_07"]]) |> sort()
		   	          ),
		   	          limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
		   	          labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
		   	        ) +
		   	        scale_fill_manual(
		   	          values = c(
		   	            "TRUE" = "black",
		   	            "FALSE" = "transparent"
		   	          )
		   	        ) +
		   	        labs(title = sample_id) )
			
			
			
		   			library(escheR)
		   			library(SpatialExperiment)
		   			library(tidyverse)
		   			library(limma)
		   			library(sessioninfo)
		   			library(here)
		   			library(spatialLIBD)
			