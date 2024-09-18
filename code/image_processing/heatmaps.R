setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')

library("here")
library("tidyr")
library("spatialLIBD")
library("dplyr")
library("ggplot2")
library("ggridges")

library(ggplot2)
library(dplyr)
library(reshape2)



readRDS(here("processed-data", "image_processing", "EDAspe.rds"))

genes_of_interest <- c("PVALB", "CLDN5", "SNAP25", "CAMK2A")
# Get the matching row indices from the rowData(spe) based on the gene_name column
matched_genes <- rownames(spe)[rowData(spe)$gene_name %in% genes_of_interest]

# Convert the dgCMatrix to a regular matrix
gene_expression_matrix <- as.matrix(assay(spe)[matched_genes, ])

# Transpose the matrix and convert it to a data frame
gene_expression_df <- as.data.frame(t(gene_expression_matrix))

# Name the columns after the genes of interest
colnames(gene_expression_df) <- genes_of_interest

# Attach the gene expression data as new columns to colData(spe)
colData(spe) <- cbind(colData(spe), gene_expression_df)

df = as.data.frame(colData(spe))


df <- df %>% mutate(PNN = ifelse(PWFA > 0.05, "PNN+", "PNN-"))
df <- df %>% mutate(DAPI = ifelse(PDAPI > 0.05 & PDAPI < 0.5 , "DAPI+", "DAPI-"))
df <- df %>% mutate(NeuN = ifelse(PNeuN > 0.05 & PNeuN < 0.3 , "NeuN+", "NeuN-"))
df <- df %>% mutate(Claudin = ifelse(PClaudin5 > 0.05 & PClaudin5 < 0.20, "Claudin+", "Claudin-"))
  
  
  # Calculate the proportion of +cells by PRECAST_07
  prop_df <- df %>%
	 group_by(PRECAST_07, sample_id, dx) %>%
	 summarize(proportion_PNN = mean(PNN == "PNN+"),
	 proportion_DAPI = mean(DAPI == "DAPI-"),
	 proportion_NeuN = mean(NeuN == "NeuN+"),
	 proportion_Claudin = mean(Claudin == "Claudin+"))
	  
  # Prioritize ntc over scz in ordering
  prop_df$dx <- factor(prop_df$dx, levels = c("ntc", "scz"))

  # Create a custom ordering for sample_id based on dx
 
  # Create the heat map
  ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_PNN)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
    labs(title = "Heat Map of PWFA by Sample ID and PRECAST_07",
         x = "Sample ID",
         y = "PRECAST_07",
         fill = "PWFA") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "right")


  # Create the heat map
   ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_NeuN)) +
     geom_tile() +
     scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
     labs(title = "Heat Map of PNeuN by Sample ID and PRECAST_07",	
          x = "Sample ID",
          y = "PRECAST_07",
          fill = "PNeuN") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1),
		            legend.position = "right")
					
					
  # Create the heat map
  ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_Claudin)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
    labs(title = "Heat Map of PClaudin5 by Sample ID and PRECAST_07",	
         x = "Sample ID",
         y = "PRECAST_07",
         fill = "PClaudin5") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1),
				  		            legend.position = "right")
	
# DAPI
 ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_DAPI)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of PDAPI- by Sample ID and PRECAST_07",	
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "PDAPI") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")

## genes

 df <- df %>% mutate(cla = ifelse(CLDN5 > 2, "cla+", "cla-"))
 df <- df %>% mutate(neuropil = ifelse(CAMK2A > 0, "neuropil+", "neuropil-"))
 df <- df %>% mutate(neuron = ifelse(SNAP25 > 2, "neuron+", "neuron-"))
 df <- df %>% mutate(pnn = ifelse(PVALB > 2, "pnn+", "pnn-"))
 
prop_df_gene <- df %>%
  group_by(PRECAST_07, sample_id, dx) %>%
  summarize(proportion_cla = mean(cla == "cla+"),
  proportion_neuropil = mean(neuropil == "neuropil+"),
  proportion_neuron = mean(neuron == "neuron+"),
  proportion_pnn = mean(pnn == "pnn+"))
 
 
  # neuropil
   ggplot(prop_df_gene, aes(x = sample_id, y = PRECAST_07, fill = proportion_neuropil)) +
     geom_tile() +
     scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
     labs(title = "Heat Map of CAMK2A by Sample ID and PRECAST_07",	
          x = "Sample ID",
          y = "PRECAST_07",
          fill = "CAMK2A") +
  	 theme_minimal() +
  	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")
 
  	
# Create the bar plot
ggplot(prop_df, aes(x = PRECAST_07, y = proportion_PNN, fill = dx)) +
  geom_bar(stat = "identity",  position = "dodge") +
  labs(title = "Proportion of PNN+ Cells by PRECAST_07",
       x = "PRECAST_07",
       y = "Proportion of PNN+ Cells") +
  theme_minimal()
	