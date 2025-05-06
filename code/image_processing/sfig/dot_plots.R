setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)

neuropil <- read_excel(here("raw-data/images/SPG_Spot_Valid/Neuropil/Supple_Table_10_Neuropil.xlsx"), sheet = "10_human_synase_markers", skip = 1)
# Remove the first 8 columns
neuropil <- neuropil[, -(1:8)]
neuropil = as.data.frame(neuropil)

genes <- c("CAMK2A",  "MAP2", "DLG4", "RBFOX3", "SNAP25", "PVALB", "CLDN5")
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% genes]
gene_expression_matrix <- as.matrix(assay(spe_ntc, "logcounts")[matched_genes, ])
gene_expression_df <- as.data.frame(t(gene_expression_matrix))
colnames(gene_expression_df) <- genes

#check
all(rownames(gene_expression_df) == rownames(colData(spe_ntc)))
gene_expression_df <- gene_expression_df[rownames(colData(spe_ntc)), ]
colData(spe_ntc) <- cbind(colData(spe_ntc), gene_expression_df)

##### barplots
df = as.data.frame(colData(spe_ntc))

colls = c("neuropil_pos", "neuropil_pos", "neuropil_pos", "neun_pos",  "neun_pos", "pnn_pos", "vasc_pos")
#colls = c("PNN", "Claudin", "NeuN", "DAPI", "DAPI", "NeuN", "DAPI")
# Loop through each gene and create the plot
pdf(here("plots", "image_processing", "enrichment", "bar.pdf"), width = 10, height = 8)

# Loop through each column in colls for x and each gene for y
for (i in c(1:7)) {
  p <- ggplot(df, aes_string(x = colls[i], y = genes[i])) +
    geom_bar(stat="identity") +
    labs(
      title = paste0(genes[i], " expression by ", colls[i]),
      x = colls[i],
      y = paste0("Expression of ", genes[i])
    ) +
    theme_minimal()
  print(p)
}

dev.off()
  

## threshold neuropil
df = as.data.frame(colData(spe_ntc))
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% neuropil[,1]]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]

df <- df %>% mutate(neuropil = ifelse(sum_expression > 1000, "neuropil+", "neuropil-"))
prop_df_gene <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_neuropil = mean(neuropil == "neuropil+"))
prop_df_gene = prop_df_gene  %>%
 mutate(
      PRECAST_07 = factor(PRECAST_07, 
		  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
           labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
           )
			
			  
prop_df <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_DAPI = mean(neuropil_pos == TRUE))
prop_df = prop_df  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          )
 
pdf(here("plots", "image_processing", "enrichment", "DAPIheatmap.pdf"), width = 10, height = 8) 
  # neuropil
p = ggplot(prop_df_gene, aes(x = sample_id, y = PRECAST_07, fill = proportion_neuropil)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of sum of expression of markergenes by Sample ID and PRECAST_07",	
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "sumexpression>1000") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")
	 print(p)
	 	
 
p = ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_DAPI)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of PDAPI by Sample ID and PRECAST_07",
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "PDAPI<0.05") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         legend.position = "right")
		 
	print(p)
	
dev.off()	 
		 
		 
## threshold neun
df = as.data.frame(colData(spe_ntc))
neun1 = read_excel(here("raw-data/images/SPG_Spot_Valid/NeuN/TableS13_marker_stats_supp_table.xlsx"), sheet="marker_stats_supp_table")
neun1 = as.data.frame(neun1)
neun1$padj = as.numeric(neun1$log.p.value)
celltypes = c("Inhib", "Excit", "Excit_L3/4/5", "Excit_L3", "Excit_L4", "Excit_L6", "Excit_L5/6", "Excit_L5", "Excit_L2/3")
neun = neun1[neun1$cellType.target %in% celltypes,]

matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_id %in% neun$gene]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]

df <- df %>% mutate(neuron = ifelse(sum_expression > 250, "neun+", "neun-"))
prop_df_gene <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_neuron = mean(neuron == "neun+"))
prop_df_gene = prop_df_gene  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          )
		  
prop_df <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_NeuN = mean(neun_pos == TRUE))
prop_df = prop_df  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          )
 
pdf(here("plots", "image_processing", "enrichment", "NeuNheatmap.pdf"), width = 10, height = 8) 
p = ggplot(prop_df_gene, aes(x = sample_id, y = PRECAST_07, fill = proportion_neuron)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of sum of expression of markergenes by Sample ID and PRECAST_07",	
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "sumexpression>250") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")
	 print(p)
	 	
 
p = ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_NeuN)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of PNeuN by Sample ID and PRECAST_07",
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "PNeuN>0.05") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         legend.position = "right")
		 
	print(p)
	
dev.off()	 

## threshold claudin
df = as.data.frame(colData(spe_ntc))
vasc = read_excel(here("raw-data/images/SPG_Spot_Valid/Vasculature/Supple_Table_2_Vas.xlsx"), sheet="Post Mortem Vascular Subcluster")
vasc = as.data.frame(vasc)
vasc$padj = as.numeric(vasc$p_val_adj)

matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% vasc$Column1]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]

df <- df %>% mutate(vasc = ifelse(sum_expression > 2500, "vasc+", "vasc-"))
prop_df_gene <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_vasc = mean(vasc == "vasc+"))
prop_df_gene = prop_df_gene  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          )
		  
prop_df <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_Claudin = mean(vasc_pos == TRUE))
prop_df = prop_df  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          )
 
pdf(here("plots", "image_processing", "enrichment", "Claudinheatmap.pdf"), width = 10, height = 8) 
p = ggplot(prop_df_gene, aes(x = sample_id, y = PRECAST_07, fill = proportion_vasc)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of sum of expression of markergenes by Sample ID and PRECAST_07",	
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "sumexpression>2500") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")
	 print(p)
	 	
 
p = ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_Claudin)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of PDAPI by Sample ID and PRECAST_07",
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "PClaudin5>0.05") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         legend.position = "right")
		 
	print(p)
	
dev.off()	

#### threshold pnn
df = as.data.frame(colData(spe_ntc))
pnn = read_excel(here("raw-data/images/SPG_Spot_Valid/PNN/Supple_Table_DataS4_PNN.xlsx"), sheet = "PNN Energy")
pnn = as.data.frame(pnn)
 
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% toupper(pnn$gene_acronym)]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]

df <- df %>% mutate(pnn = ifelse(sum_expression > 5000, "pnn+", "pnn-"))
prop_df_gene <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_pnn = mean(pnn == "pnn+"))
prop_df_gene = prop_df_gene  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          )
		  
prop_df <- df %>%
group_by(PRECAST_07, sample_id) %>%
summarize(proportion_WFA = mean(pnn_pos == TRUE))
prop_df = prop_df  %>%
mutate(
     PRECAST_07 = factor(PRECAST_07, 
	  levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
          labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
          ) 
 
pdf(here("plots", "image_processing", "enrichment", "WFAheatmap.pdf"), width = 10, height = 8) 
p = ggplot(prop_df_gene, aes(x = sample_id, y = PRECAST_07, fill = proportion_pnn)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of sum of expression of markergenes by Sample ID and PRECAST_07",	
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "sumexpression>5000") +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")
	 print(p)
	 	
 
p = ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = proportion_WFA)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Heat Map of PDAPI by Sample ID and PRECAST_07",
        x = "Sample ID",
        y = "PRECAST_07",
        fill = "PWFA>0.05") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         legend.position = "right")
		 
	print(p)
	
dev.off()	



##### dot plots ######
spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
df = as.data.frame(colData(spe_ntc))

neuropil <- read_excel(here("raw-data/images/SPG_Spot_Valid/Neuropil/Supple_Table_10_Neuropil.xlsx"), sheet = "10_human_synase_markers", skip = 1)
neuropil <- neuropil[, -(1:8)]
neuropil = as.data.frame(neuropil)
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% neuropil[,1]]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression_neuropil <- sum_expression[rownames(df)]

neun1 = read_excel(here("raw-data/images/SPG_Spot_Valid/NeuN/TableS13_marker_stats_supp_table.xlsx"), sheet="marker_stats_supp_table")
neun1 = as.data.frame(neun1)
neun1$padj = as.numeric(neun1$log.p.value)
celltypes = c("Inhib", "Excit", "Excit_L3/4/5", "Excit_L3", "Excit_L4", "Excit_L6", "Excit_L5/6", "Excit_L5", "Excit_L2/3")
neun = neun1[neun1$cellType.target %in% celltypes,]
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_id %in% neun$gene]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression_neun <- sum_expression[rownames(df)]

vasc = read_excel(here("raw-data/images/SPG_Spot_Valid/Vasculature/Supple_Table_2_Vas.xlsx"), sheet="Post Mortem Vascular Subcluster")
vasc = as.data.frame(vasc)
vasc$padj = as.numeric(vasc$p_val_adj)
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% vasc$Column1]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression_vasc <- sum_expression[rownames(df)]

pnn = read_excel(here("raw-data/images/SPG_Spot_Valid/PNN/Supple_Table_DataS4_PNN.xlsx"), sheet = "PNN Energy")
pnn = as.data.frame(pnn)
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% toupper(pnn$gene_acronym)]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression_pnn <- sum_expression[rownames(df)]

df <- df %>% mutate(vasc = ifelse(sum_expression_vasc > 2000, TRUE, FALSE))
df <- df %>% mutate(neun = ifelse(sum_expression_neun > 500, TRUE, FALSE))
df <- df %>% mutate(neuropil = ifelse(sum_expression_neuropil <2500, TRUE, FALSE))
df <- df %>% mutate(pnn = ifelse(sum_expression_pnn > 2000, TRUE, FALSE))

library(ggplot2)

# Create a PDF file to save all plots
pdf(here("plots", "image_processing", "enrichment", "sum_expression_plots.pdf"), width = 8, height = 6)

# Define the pairs of columns for plotting
expression_columns <- c("sum_expression_neuropil", "sum_expression_vasc", "sum_expression_pnn", "sum_expression_neun")
pos_columns <- c("neuropil_pos", "vasc_pos", "pnn_pos", "neun_pos")

# Loop through each pair of columns
for (i in seq_along(expression_columns)) {
  sum_expression_col <- expression_columns[i]
  pos_col <- pos_columns[i]
  
  # Dynamically create the plot
  p <- ggplot(df, aes_string(x = "sample_id", y = sum_expression_col, color = pos_col)) +
    geom_jitter(size = 3, alpha = 0.5, width = 0.2, height = 0) +
    labs(
      title = paste("Scatter Plot of", sum_expression_col),
      x = "Sample ID",
      y = sum_expression_col
    ) +
    theme_minimal() +
    scale_color_viridis_d(option = "plasma")
  
  # Print the plot to the PDF file
  print(p)
}

dev.off()

prop_df <- df %>%
group_by(PRECAST_07) %>%
summarize(proportion_pnn = mean(pnn_pos == TRUE),
	      mean_pnn = mean(sum_expression_pnn),
		  prop_pnn = mean(pnn==pnn_pos),
		  proportion_neuropil = mean(neuropil_pos == TRUE),
		  mean_neuropil = mean(sum_expression_neuropil),
		  prop_neuropil = mean(neuropil==neuropil_pos),
		  proportion_neun = mean(neun_pos == TRUE),
		  mean_neun = mean(sum_expression_neun),
		  prop_neun = mean(neun==neun_pos),
		  proportion_vasc = mean(vasc_pos == TRUE),
		  mean_vasc = mean(sum_expression_vasc),
		  prop_vasc = mean(vasc==vasc_pos),
) %>%
mutate(
    PRECAST_07 = factor(PRECAST_07, 
                        levels = c("spd07", "spd06", "spd02", "spd05", "spd03", "spd01", "spd04"),
                        labels = c("spd07-L1", "spd06-L2/3", "spd02-L3/4", "spd05-L5", "spd03-L6", "spd01-L6/WM", "spd04-WM"))
  )


properties <- c("pnn", "neuropil", "neun", "vasc")
pdf(here("plots", "image_processing", "enrichment", "dot_plots.pdf"), width = 8, height = 6)

for (prop in properties) {
  p = ggplot(prop_df, aes_string(
    x = "PRECAST_07", 
    y = paste0("proportion_", prop), 
    fill = paste0("mean_", prop), 
    size = paste0("proportion_", prop)
  )) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_fill_gradient(low = "white", high = "black") +
    labs(
      title = paste("Dot Plot of", toupper(prop)),
      x = "PRECAST_07",
      y = paste("Proportion of matching spots"),
      fill = "Mean Expression",
      size = paste("Proportion", toupper(prop) , "+spots")
    ) +
    theme_minimal()
    print(p) # Print plot to the PDF
}

dev.off()

unique_spd_values <- unique(df$PRECAST_07)

# Initialize empty lists to store results
prop_neuropil <- list()
prop_neun <- list()
prop_pnn <- list()
prop_vasc <- list()

# Loop through each unique PRECAST_07 value and calculate the proportion
for (spd in unique_spd_values) {
  # Calculate the proportion for neuropil
  proportion_neuropil <- sum(df$neuropil_pos == TRUE & df$PRECAST_07 == spd) / sum(df$neuropil_pos == TRUE)
  prop_neuropil[[spd]] <- proportion_neuropil
  mean_neuropil[[spd]] = mean(df$sum_expression_neuropil[df$PRECAST_07 == spd])
  # Calculate the proportion for NeuN
  proportion_neun <- sum(df$neun_pos == TRUE & df$PRECAST_07 == spd) / sum(df$neun_pos == TRUE)
  prop_neun[[spd]] <- proportion_neun
  mean_neun[[spd]] = mean(df$sum_expression_neun[df$PRECAST_07 == spd])
  # Calculate the proportion for PNN
  proportion_pnn <- sum(df$pnn_pos == TRUE & df$PRECAST_07 == spd) / sum(df$pnn_pos == TRUE)
  prop_pnn[[spd]] <- proportion_pnn
  mean_pnn[[spd]] = mean(df$sum_expression_pnn[df$PRECAST_07 == spd])
  # Calculate the proportion for vascular position
  proportion_vasc <- sum(df$vasc_pos == TRUE & df$PRECAST_07 == spd) / sum(df$vasc_pos == TRUE)
  prop_vasc[[spd]] <- proportion_vasc
  mean_vasc[[spd]] = mean(df$sum_expression_vasc[df$PRECAST_07 == spd])
}

# Combine the results into a data frame
prop_df1 <- data.frame(
  PRECAST_07 = unique_spd_values,
  Neuropil = unlist(prop_neuropil),
  NeuropilM = unlist(mean_neuropil),
  NeuN = unlist(prop_neun),
  NeuNM = unlist(mean_neun),
  WFA = unlist(prop_pnn),
  WFAM = unlist(mean_pnn),
  Claudin_5 = unlist(prop_vasc),
  Claudin_5M = unlist(mean_vasc)
)

properties <- c("WFA", "Neuropil", "NeuN", "Claudin_5")
pdf(here("plots", "image_processing", "enrichment", "dot_plots1.pdf"), width = 8, height = 6)

for (prop in properties) {
  p = ggplot(prop_df1, aes_string(
    x = "PRECAST_07", 
    y =  prop, 
    fill = paste0("mean_", prop), 
    size = paste0("proportion_", prop)
  )) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_fill_gradient(low = "white", high = "black") +
    labs(
      title = paste("Dot Plot of", toupper(prop)),
      x = "PRECAST_07",
      y = paste("Proportion of matching spots"),
      fill = "Mean Expression",
      size = paste("Proportion", toupper(prop) , "+spots")
    ) +
    theme_minimal()
    print(p) # Print plot to the PDF
}

dev.off()