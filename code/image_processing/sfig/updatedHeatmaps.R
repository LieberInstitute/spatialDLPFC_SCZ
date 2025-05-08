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
df = as.data.frame(colData(spe_ntc))

neuropil <- read_excel(here("raw-data/images/SPG_Spot_Valid/Neuropil/Supple_Table_10_Neuropil.xlsx"), sheet = "10_human_synase_markers", skip = 1)
# Remove the first 8 columns
neuropil <- neuropil[, -(1:8)]
neuropil = as.data.frame(neuropil)
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% neuropil[,1]]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]
df <- df %>% mutate(neuropil = ifelse(sum_expression > 1000, TRUE, FALSE))

neun1 = read_excel(here("raw-data/images/SPG_Spot_Valid/NeuN/TableS13_marker_stats_supp_table.xlsx"), sheet="marker_stats_supp_table")
neun1 = as.data.frame(neun1)
neun1$padj = as.numeric(neun1$log.p.value)
celltypes = c("Inhib", "Excit", "Excit_L3/4/5", "Excit_L3", "Excit_L4", "Excit_L6", "Excit_L5/6", "Excit_L5", "Excit_L2/3")
neun = neun1[neun1$cellType.target %in% celltypes,]
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_id %in% neun$gene]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]
df <- df %>% mutate(neuron = ifelse(sum_expression > 250, TRUE, FALSE))

vasc = read_excel(here("raw-data/images/SPG_Spot_Valid/Vasculature/Supple_Table_2_Vas.xlsx"), sheet="Post Mortem Vascular Subcluster")
vasc = as.data.frame(vasc)
vasc$padj = as.numeric(vasc$p_val_adj)
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% vasc$Column1]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]
df <- df %>% mutate(vasc = ifelse(sum_expression > 2500, TRUE, FALSE))

pnn = read_excel(here("raw-data/images/SPG_Spot_Valid/PNN/Supple_Table_DataS4_PNN.xlsx"), sheet = "PNN Energy")
pnn = as.data.frame(pnn)
matched_genes <- rownames(spe_ntc)[rowData(spe_ntc)$gene_name %in% toupper(pnn$gene_acronym)]
expr_data <- assay(spe_ntc, "counts")[matched_genes, ]
sum_expression <- colSums(expr_data)
df$sum_expression <- sum_expression[rownames(df)]
df <- df %>% mutate(pnn = ifelse(sum_expression > 5000, TRUE, FALSE))

# Get unique values
unique_samples <- unique(df$sample_id)
unique_spd_values <- unique(df$PRECAST_07)

# Initialize empty lists
prop_neuropil <- list()
prop_neun <- list()
prop_pnn <- list()
prop_vasc <- list()

# Loop through each sample and each PRECAST_07 domain
for (samp in unique_samples) {
  for (spd in unique_spd_values) {
    subset_df <- df[df$sample_id == samp & df$PRECAST_07 == spd, ]
    
    # Store key name
    key <- paste(samp, spd, sep = "_")
    
    # Proportion and mean for neuropil
    prop_neuropil[[key]] <- sum(subset_df$neuropil, na.rm = TRUE) / sum(df$neuropil == TRUE & df$sample_id == samp, na.rm = TRUE)
    
    # Proportion and mean for NeuN
    prop_neun[[key]] <- sum(subset_df$neuron, na.rm = TRUE) / sum(df$neuron == TRUE & df$sample_id == samp, na.rm = TRUE)
    
    # Proportion and mean for PNN
    prop_pnn[[key]] <- sum(subset_df$pnn, na.rm = TRUE) / sum(df$pnn == TRUE & df$sample_id == samp, na.rm = TRUE)
    
    # Proportion and mean for vascular
    prop_vasc[[key]] <- sum(subset_df$vasc, na.rm = TRUE) / sum(df$vasc == TRUE & df$sample_id == samp, na.rm = TRUE)
  }
}

keys <- names(prop_neuropil)

# Use regular expressions to split on the last underscore
sample_id <- sub("_(spd[0-9]+)$", "", keys)
PRECAST_07 <- sub("^.*_(spd[0-9]+)$", "\\1", keys)

# Combine into a data frame
gene_df <- data.frame(
  sample_id = sample_id,
  PRECAST_07 = PRECAST_07,
  P_neuropil = unlist(prop_neuropil),
  P_neun = unlist(prop_neun),
  P_pnn = unlist(prop_pnn),
  P_vasc = unlist(prop_vasc)
)%>%
  mutate(PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
     labels = c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3","SpD07-L1" ))
  )
  
# Initialize empty lists
prop_neuropil <- list()
prop_neun <- list()
prop_pnn <- list()
prop_vasc <- list()

# Loop through each sample and each PRECAST_07 domain
for (samp in unique_samples) {
  for (spd in unique_spd_values) {
    subset_df <- df[df$sample_id == samp & df$PRECAST_07 == spd, ]
    
    # Store key name
    key <- paste(samp, spd, sep = "_")
    
    # Proportion and mean for neuropil
    prop_neuropil[[key]] <- sum(subset_df$neuropil_pos, na.rm = TRUE) / sum(df$neuropil_pos == TRUE & df$sample_id == samp, na.rm = TRUE)
    
    # Proportion and mean for NeuN
    prop_neun[[key]] <- sum(subset_df$neun_pos, na.rm = TRUE) / sum(df$neun_pos == TRUE & df$sample_id == samp, na.rm = TRUE)
    
    # Proportion and mean for PNN
    prop_pnn[[key]] <- sum(subset_df$pnn_pos, na.rm = TRUE) / sum(df$pnn_pos == TRUE & df$sample_id == samp, na.rm = TRUE)
    
    # Proportion and mean for vascular
    prop_vasc[[key]] <- sum(subset_df$vasc_pos, na.rm = TRUE) / sum(df$vasc_pos == TRUE & df$sample_id == samp, na.rm = TRUE)
  }
}

keys <- names(prop_neuropil)

# Use regular expressions to split on the last underscore
sample_id <- sub("_(spd[0-9]+)$", "", keys)
PRECAST_07 <- sub("^.*_(spd[0-9]+)$", "\\1", keys)

# Combine into a data frame
prop_df <- data.frame(
  sample_id = sample_id,
  PRECAST_07 = PRECAST_07,
  P_neuropil = unlist(prop_neuropil),
  P_neun = unlist(prop_neun),
  P_pnn = unlist(prop_pnn),
  P_vasc = unlist(prop_vasc)
)%>%
  mutate(PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
     labels = c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3","SpD07-L1" ))
  )
  
  
pdf(here("plots", "image_processing", "enrichment", "vascheatmap_update.pdf"), width = 10, height = 8) 
  # neuropil
p = ggplot(gene_df, aes(x = sample_id, y = PRECAST_07, fill = P_vasc)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Proportion of Vasculature+ spots with sum of marker gene expression > 2500",	
        x = NULL,
        y = NULL,
		fill = NULL) +
	 theme_minimal() +
	 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 20), legend.position = "right", 
	 axis.text.y = element_text(size = 20), plot.title = element_text(hjust = 0.8, size = 20))

	 print(p)
	 
 
p = ggplot(prop_df, aes(x = sample_id, y = PRECAST_07, fill = P_vasc)) +
   geom_tile() +
   scale_fill_gradient(low = "white", high = "black", na.value = "grey50") +
   labs(title = "Proportion of Claudin5+ spots with segmented Claudin5 signal > 5% of spot area",
        x = NULL,
        y = NULL,
		fill = NULL) +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 20), legend.position = "right", 
   axis.text.y = element_text(size = 20), plot.title = element_text(hjust = 0.7, size = 20))
		 
	print(p)
	
dev.off()	 
		 