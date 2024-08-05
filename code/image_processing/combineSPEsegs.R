setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
library("readxl")
library("here")
library("tidyr")
library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("ggridges")
#spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_wo_spg_N63.rds"))
#PRECAST_df <- readRDS(here("processed-data/rds/spatial_cluster","PRECAST", "test_clus_label_df_semi_inform_k_2-16.rds"))
#stopifnot(nrow(PRECAST_df) == ncol(spe))
#
#precast_vars <- grep("^PRECAST_", colnames(PRECAST_df),value = TRUE)
#spe <- spe[, spe$key %in% PRECAST_df$key]
## raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
#col_data_df <- PRECAST_df |> right_join(
#    colData(spe) |> data.frame(),
#    by = c("key"),
#    relationship = "one-to-one"
#  )
#
#rownames(col_data_df) <- colnames(spe)
#colData(spe) <- DataFrame(col_data_df)


spe = readRDS(here("processed-data/rds/spatial_cluster/PRECAST/spe_wo_spg_N63_PRECAST.rds"))
samples = read_excel(here("raw-data","experiment_info","VisiumSPG_PNN_Master.xlsx"))
samples <- as.data.frame(samples %>% select('Slide #', 'Array #'))
samples$sample_id = paste0(samples[,1],"_",samples[,2])
samples$SPpath = paste0('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/',samples$sample_id,'/outs/spatial/tissue_spot_counts.csv')
segmentations_list <-
  lapply(samples$sample_id, function(sample_id) {
    file <-samples$SPpath[samples$sample_id == sample_id]
    if (!file.exists(file)) {
      return(NULL)
    }
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sample_id)
    return(x)
  })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
  Reduce(function(...) {
    merge(..., all = TRUE)
  }, segmentations_list[lengths(segmentations_list) > 0])


## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spe) <- cbind(colData(spe), segmentation_info)

saveRDS(here("processed-data", "image_processing", "EDAspe.rds"), spe)


coldata_df = as.data.frame(colData(spe))

p1 = ggplot(coldata_df, aes(x = PWFA)) + geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
    labs(title = "Distribution of WFA", x = "PWFA", y = "Frequency") +
    theme_minimal()
	  
p11 = ggplot(coldata_df, aes(x = PWFA)) + geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
	    labs(title = "Distribution of WFA", x = "PWFA >= 2%", y = "Frequency") +
	    theme_minimal()+xlim(0.001,1)+geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", size = 1)
		
p12 = ggplot(coldata_df, aes(y = PWFA, x = sample_id, fill = dx)) +geom_violin(trim = FALSE, alpha = 0.7) +
		labs(title = "Distribution of WFA", x = "PWFA >=5%", y = "Frequency") +
		theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +ylim(0.05,1) +facet_wrap(~ dx, scales = "free")

p2 = ggplot(coldata_df, aes(x = dx, y = PWFA)) + geom_boxplot(outlier.shape = NA) +
	labs(title = "Distribution of WFA by Dx",x = "Dx",y = "PWFA") + geom_jitter(width = 0.2, alpha = 0.5) +
	theme_minimal()

p21 = ggplot(coldata_df, aes(x = dx, y = PWFA)) + geom_boxplot(outlier.shape = NA) +
		labs(title = "Distribution of WFA by Dx",x = "Dx",y = "PWFA >= 2%") + geom_jitter(width = 0.2, alpha = 0.5) +
		theme_minimal()+ylim(0.02,1)

p22 = ggplot(coldata_df, aes(x = dx, y = PWFA, fill = dx)) + geom_violin(trim = FALSE, alpha = 0.7) +
		labs(title = "Distribution of WFA by Dx",x = "Dx",y = "PWFA >= 25%") + 
		theme_minimal()+ylim(0.25,1) + theme(legend.position = "none") 

	  	
p3 = ggplot(coldata_df, aes(x = sex, y = PWFA)) + geom_boxplot(outlier.shape = NA) +
	labs(title = "Distribution of WFA by sex",x = "sex",y = "PWFA") + geom_jitter(width = 0.2, alpha = 0.5) +
	theme_minimal()

p31 = ggplot(coldata_df, aes(x = sex, y = PWFA)) + geom_boxplot(outlier.shape = NA) +
		labs(title = "Distribution of WFA by sex",x = "sex",y = "PWFA >= 25%") + geom_jitter(width = 0.2, alpha = 0.5) +
		theme_minimal()+ylim(0.25,1)

p32 = ggplot(coldata_df, aes(x = sex, y = PWFA, fill = dx)) + geom_violin(trim = FALSE, alpha = 0.7) +
		labs(title = "Distribution of WFA by sex",x = "sex",y = "PWFA >= 2%") + 
				theme_minimal()+ylim(0.02,1) 

			    pdf(here("plots","image_processing","pnn_EDA_PWFA_M.pdf"), width = 10, height = 8)
			    print(p1)
				print(p11)
				print(p12)
				print(p2)
				print(p21)
				print(p22)
				print(p3)
				print(p31)
				print(p32)
			    dev.off()
					
coldata_df$IWFA <- ifelse(is.na(coldata_df$IWFA), 0, coldata_df$IWFA)
filtered_data <- coldata_df[coldata_df$IWFA < 1, ]

p1 = ggplot(filtered_data, aes(x = IWFA)) + geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
    labs(title = "Distribution of WFA", x = "IWFA", y = "Frequency") +
    theme_minimal()
	  
p11 = ggplot(filtered_data, aes(x = IWFA)) + geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
	    labs(title = "Distribution of WFA", x = "IWFA >= 0", y = "Frequency") +
	    theme_minimal()+xlim(0.01,1)+geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", size = 1)

p12 = ggplot(coldata_df, aes(y = IWFA, x = sample_id, fill = dx)) +geom_violin(trim = FALSE, alpha = 0.7) +
	labs(title = "Distribution of WFA", x = "IWFA>0", y = "Frequency") +
	 theme_minimal()+ylim(0.1,3)+facet_wrap(~ dx, scales = "free")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
			
p2 = ggplot(filtered_data, aes(x = dx, y = IWFA)) + geom_boxplot(outlier.shape = NA) +
	labs(title = "Distribution of WFA by Dx",x = "Dx",y = "IWFA") + geom_jitter(width = 0.2, alpha = 0.5) +
	theme_minimal()

p21 = ggplot(filtered_data, aes(x = dx, y = IWFA)) + geom_boxplot(outlier.shape = NA) +
		labs(title = "Distribution of WFA by Dx",x = "Dx",y = "IWFA >0") + geom_jitter(width = 0.2, alpha = 0.5) +
		theme_minimal()+ylim(0.1,1)

p22 = ggplot(filtered_data, aes(x = dx, y = IWFA, fill = dx)) + geom_violin(trim = FALSE, alpha = 0.7) +
		labs(title = "Distribution of WFA by Dx",x = "Dx",y = "IWFA >0") + 
		theme_minimal()+ylim(0.1,1) + theme(legend.position = "none") 
		
	  	
p3 = ggplot(filtered_data, aes(x = sex, y = IWFA)) + geom_boxplot(outlier.shape = NA) +
	labs(title = "Distribution of WFA by sex",x = "sex",y = "IWFA") + geom_jitter(width = 0.2, alpha = 0.5) +
	theme_minimal()

p31 = ggplot(filtered_data, aes(x = sex, y = IWFA)) + geom_boxplot(outlier.shape = NA) +
		labs(title = "Distribution of WFA by sex",x = "sex",y = "IWFA >0") + geom_jitter(width = 0.2, alpha = 0.5) +
		theme_minimal()+ylim(0.1,1)

p32 = ggplot(filtered_data, aes(x = sex, y = IWFA, fill = dx)) + geom_violin(trim = FALSE, alpha = 0.7) +
		labs(title = "Distribution of WFA by sex",x = "sex",y = "IWFA >0") + 
				theme_minimal()+ylim(0.1,1) 
	
			    pdf(here("plots","image_processing","pnn_EDA_IWFA_M.pdf"), width = 10, height = 8)
			    print(p1)
				print(p11)
				print(p12)
				print(p2)
				print(p21)
				print(p22)
				print(p3)
				print(p31)
				print(p32)
			    dev.off()

gene <- rownames(spe)[match("PVALB", rowData(spe)$gene_name)]		  
spT <- spe[gene, ]
spT$gene = as.numeric(t(assays(spT)$counts))
coldata_df = as.data.frame(colData(spT))
coldata_df_ntc = coldata_df[coldata_df$dx == 'ntc', ]

proportion_data <- coldata_df_ntc %>%
  group_by(sample_id, PRECAST_07, CNWFA_GT_0 = PWFA > 0.05) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(sample_id, PRECAST_07) %>%
  mutate(Proportion = Count / sum(Count))

	# Create stacked bar plot
p1 = ggplot(proportion_data, aes(x = PRECAST_07, y = Proportion, fill = CNWFA_GT_0)) +
	  geom_bar(stat = "identity", position = "stack") +
	  labs(title = "Proportion of WFA+ by spD in NTCs", x = "Sample ID", y = "Proportion", fill = "PWFA > 0.05") +
	  theme_minimal()+facet_wrap(~ sample_id)+
	  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

proportion_data <- coldata_df_ntc %>%
  group_by(sample_id, PRECAST_07, PVALB = gene > 1) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(sample_id, PRECAST_07) %>%
  mutate(Proportion = Count / sum(Count))

p2 = ggplot(proportion_data, aes(x = PRECAST_07, y = Proportion, fill = PVALB)) +
          geom_bar(stat = "identity", position = "stack") +
          labs(title = "Proportion of PVALB+ by spD in NTCs", x = "Sample ID", y = "Proportion", fill = "PVALB >   1") +
          theme_minimal()+facet_wrap(~ sample_id)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

		  pdf(here("plots","image_processing","pnnVSpvalb_ntc_M.pdf"), width = 10, height = 8)
		  print(p1)
		  print(p2)
		  dev.off()	  	  

coldata_df_scz = coldata_df[coldata_df$dx == 'scz', ]		  
proportion_data <- coldata_df_scz %>%
	    group_by(sample_id, PRECAST_07, CNWFA_GT_0 = PWFA > 0.05) %>%
	    summarise(Count = n()) %>%
	    ungroup() %>%
	    group_by(sample_id, PRECAST_07) %>%
	    mutate(Proportion = Count / sum(Count))

	  	# Create stacked bar plot
p1 = ggplot(proportion_data, aes(x = PRECAST_07, y = Proportion, fill = CNWFA_GT_0)) +
	  geom_bar(stat = "identity", position = "stack") +
	  labs(title = "Proportion of WFA+ by spD in SCZs", x = "Sample ID", y = "Proportion", fill = "PWFA > 0.05") +
	  theme_minimal()+facet_wrap(~ sample_id)+
	  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
	  
	  		  
proportion_data <- coldata_df_scz %>%
  group_by(sample_id, PRECAST_07, PVALB = gene > 1) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(sample_id, PRECAST_07) %>%
  mutate(Proportion = Count / sum(Count))
  
  p2 = ggplot(proportion_data, aes(x = PRECAST_07, y = Proportion, fill = PVALB)) +
          geom_bar(stat = "identity", position = "stack") +
          labs(title = "Proportion of PVALB+ by spD in scz", x = "Sample ID", y = "Proportion", fill = "PVALB >   1") +
          theme_minimal()+facet_wrap(~ sample_id)+
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
                  
				  pdf(here("plots","image_processing","pnnVSpvalb_scz_M.pdf"), width = 10, height = 8)
				  print(p1)
				  print(p2)
				  dev.off()	  	  
library("ggpmisc")

r_squared_df <- coldata_df %>%
  group_by(PRECAST_07,dx) %>%
  summarize(corr = round(cor(PWFA, gene), 2),)
  
  r_squared_df <- coldata_df %>%
    group_by(PRECAST_07) %>%
    summarize(corr = round(cor(PWFA, gene), 2),Su = )
  
 		  
	  
	proportion_data <- coldata_df %>%
		    group_by(dx,PRECAST_07, CNWFA_GT_0 = PWFA > 0.10) %>%
		    summarise(Count = n()) %>%
		    ungroup() %>%
		    group_by(dx,PRECAST_07) %>%
		    mutate(Proportion = Count / sum(Count))

		  	# Create stacked bar plot
	p1 = ggplot(proportion_data, aes(x = PRECAST_07, y = Proportion, fill = CNWFA_GT_0)) +
		  geom_bar(stat = "identity", position = "stack") +
		  labs(title = "Proportion of WFA+ by spD", x = "PRECAST domains", y = "Proportion", fill = "PWFA > 10%") +
		  theme_minimal()+facet_wrap(~ dx)+
		  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
		  
pdf(here("plots","image_processing","pnnVSspD.pdf"), width = 10, height = 8)
print(p1)
dev.off()