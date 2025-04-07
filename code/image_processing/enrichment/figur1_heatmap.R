setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  #library(spatialLIBD)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(reshape2)
  library(gridExtra)
})

#spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_w_spg_N63.rds"))
#
### Load SpD data ----
#finalized_spd <- readRDS(here("processed-data/rds/spatial_cluster", "PRECAST","test_clus_label_df_semi_inform_k_2-16.rds"))
#
### Attach SpD label to spe ----
#col_data_df <- colData(spe) |>
#  data.frame() |>
#  left_join(
#    finalized_spd,
#    by = c("key"),
#    relationship = "one-to-one"
#  )
#
#rownames(col_data_df) <- colnames(spe)
#colData(spe) <- DataFrame(col_data_df)
#
## Call SPG spots ----
#spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
## NOTE: neuropil spot are spots doesn't have DAPI staining
#spe$neuropil_pos <- ifelse(spe$spg_PDAPI > 0.05,FALSE, TRUE)
#spe$neun_pos <- ifelse(spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,TRUE, FALSE)
#spe$vasc_pos <- ifelse(spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,TRUE, FALSE)
#
#spe_ntc <- spe[, colData(spe)$dx == "ntc"]
#spe_scz <- spe[, colData(spe)$dx == "scz"]
#
#saveRDS(spe_ntc, here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
#saveRDS(spe_scz, here("processed-data", "image_processing", "enrichment", "spe_scz.rds"))


spe_ntc = readRDS(here("processed-data/image_processing/enrichment/spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
df = as.data.frame(colData(spe_ntc))

# Create the heatmap layer
prop_df <- df %>%
group_by(PRECAST_07) %>%
summarize(Neuropil = mean(neuropil_pos == TRUE), NeuNs = mean(neun_pos == TRUE), 
WFA = mean(pnn_pos == TRUE), Claudin5 = mean(!(pnn_pos | neuropil_pos | neun_pos | vasc_pos)))

prop_long <- prop_df %>%
  pivot_longer(cols = -PRECAST_07, names_to = "Category", values_to = "Proportion") %>%
  mutate(
    Category = factor(Category, levels = c("Neuropil", "NeuN", "WFA", "Claudin5")),
    PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
                        labels = c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3","SpD07-L1" ))
  )
	
pdf(here("plots", "image_processing", "enrichment", "ntc_layer.pdf"), width = 8, height = 6)
ggplot(prop_long, aes(x = Category, y = PRECAST_07, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = "Category",
    y = "PRECAST_07",
    title = "Heatmap of Proportions"
  )
dev.off()  

# Create the heatmap by diagnosis

prop_df <- df %>%
group_by(PRECAST_07) %>%
summarize(Neuropil = sum(neuropil_pos == TRUE)/dim(df)[1], NeuN = sum(neun_pos == TRUE)/dim(df)[1], 
WFA = sum(pnn_pos == TRUE)/dim(df)[1], Claudin5 = sum(vasc_pos == TRUE)/dim(df)[1], 
NONE = sum(!(pnn_pos | neuropil_pos | neun_pos | vasc_pos))/dim(df)[1])

prop_long <- prop_df %>%
  pivot_longer(cols = -PRECAST_07, names_to = "Category", values_to = "Proportion") %>%
  mutate(
    Category = factor(Category, levels = c("Neuropil", "NeuN", "WFA", "Claudin5", "NONE")),
    PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
                        labels = c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3","SpD07-L1" ))
  )
	
pdf(here("plots", "image_processing", "enrichment", "ntc_allsample.pdf"), width = 8, height = 6)
ggplot(prop_long, aes(x = Category, y = PRECAST_07, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "Proportion") +
  geom_text(aes(label = scales::percent(round(Proportion,2))), color = "black", size = 3) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = "Category",
    y = "PRECAST_07",
    title = "Heatmap of Proportions"
  )
dev.off()  

# Create the heatmap by channel

# Create a vector of unique PRECAST_07 values
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
  
    # Calculate the proportion for NeuN
    proportion_neun <- sum(df$neun_pos == TRUE & df$PRECAST_07 == spd) / sum(df$neun_pos == TRUE)
    prop_neun[[spd]] <- proportion_neun
    
    # Calculate the proportion for PNN
    proportion_pnn <- sum(df$pnn_pos == TRUE & df$PRECAST_07 == spd) / sum(df$pnn_pos == TRUE)
    prop_pnn[[spd]] <- proportion_pnn
   
    # Calculate the proportion for vascular position
    proportion_vasc <- sum(df$vasc_pos == TRUE & df$PRECAST_07 == spd) / sum(df$vasc_pos == TRUE)
    prop_vasc[[spd]] <- proportion_vasc
 }

 # Combine the results into a data frame
 prop_df <- data.frame(
    PRECAST_07 = unique_spd_values,
    Neuropil = unlist(prop_neuropil),
    NeuN = unlist(prop_neun),
    WFA = unlist(prop_pnn),
    Claudin_5 = unlist(prop_vasc)
)


# Melt the dataframe for plotting
prop_melted <- melt(prop_df, id.vars = "PRECAST_07") %>%
  mutate(
     PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
                      labels = c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3","SpD07-L1" ))
  )

# Plot the heatmap
plot1 <- ggplot(prop_melted[prop_melted$variable == "Neuropil",], aes(x = variable, y = PRECAST_07, fill = value)) +
  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "black", limits = c(0,1))+  
  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill = guide_legend(title = NULL)) + 
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
plot2 <- ggplot(prop_melted[prop_melted$variable == "NeuN",], aes(x = variable, y = PRECAST_07, fill = value)) +
  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "yellow", limits = c(0, 1)) + 
  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill = guide_legend(title = NULL)) +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
plot3 <- ggplot(prop_melted[prop_melted$variable == "WFA",], aes(x = variable, y = PRECAST_07, fill = value)) +
  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "magenta", limits = c(0, 1)) + 
  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill = guide_legend(title = NULL)) +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
plot4 <- ggplot(prop_melted[prop_melted$variable == "Claudin_5",], aes(x = variable, y = PRECAST_07, fill = value)) +
  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) + 
  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill = guide_legend(title = NULL)) +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3)
ggsave("combined_plot.png", plot = plotC, width = 12, height = 10, dpi = 300)

pdf(here("plots", "image_processing", "enrichment", "ntc_combined.pdf"), width = 12, height = 10)
ggplot(prop_melted, aes(x = variable, y = PRECAST_07, fill = value)) +
  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "black", limits = c(0,1))+  
  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill = guide_legend(title = NULL)) + 
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
dev.off()


###################################
####### scz heatmap #########
spe_scz = readRDS(here("processed-data/image_processing/enrichment/spe_scz.rds"))
colData(spe_scz)$slide_id <- sapply(strsplit(colData(spe_scz)$sample_id, "_"), `[`, 1)
df = as.data.frame(colData(spe_scz))

prop_df <- df %>%
group_by(PRECAST_07) %>%
summarize(Neuropil_spots = mean(neuropil_pos == TRUE), NeuN_spots = mean(neun_pos == TRUE), 
WFA_spots = mean(pnn_pos == TRUE), Claudin5_spots = mean(vasc_pos == TRUE))

prop_long <- prop_df %>%
  pivot_longer(cols = -PRECAST_07, names_to = "Category", values_to = "Proportion") %>%
  mutate(
    Category = factor(Category, levels = c("Neuropil_spots", "NeuN_spots", "WFA_spots", "Claudin5_spots")),
    PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
                        labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
  )
  # Create the heatmap
  pdf(here("plots", "image_processing", "enrichment", "scz_layer.pdf"), width = 8, height = 6)
  ggplot(prop_long, aes(x = Category, y = PRECAST_07, fill = Proportion)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "black", name = "Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    labs(
      x = "Category",
      y = "PRECAST_07",
      title = "Heatmap of Proportions"
    )
	dev.off()
	
	# Create the heatmap diagnosis

	prop_df <- df %>%
	group_by(PRECAST_07) %>%
	summarize(Neuropil_spots = sum(neuropil_pos == TRUE)/dim(df)[1], NeuN_spots = sum(neun_pos == TRUE)/dim(df)[1], 
	WFA_spots = sum(pnn_pos == TRUE)/dim(df)[1], Claudin5_spots = sum(vasc_pos == TRUE)/dim(df)[1])

	prop_long <- prop_df %>%
	  pivot_longer(cols = -PRECAST_07, names_to = "Category", values_to = "Proportion") %>%
	  mutate(
	    Category = factor(Category, levels = c("Neuropil_spots", "NeuN_spots", "WFA_spots", "Claudin5_spots")),
	    PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
	                        labels = c("spd04-WM", "spd01-L6/WM", "spd03-L6", "spd05-L5", "spd02-L3/4", "spd06-L2/3","spd07-L1" ))
	  )
	
	pdf(here("plots", "image_processing", "enrichment", "scz_allsample.pdf"), width = 8, height = 6)
	ggplot(prop_long, aes(x = Category, y = PRECAST_07, fill = Proportion)) +
	  geom_tile(color = "white") +
	  scale_fill_gradient(low = "white", high = "black", name = "Proportion") +
	  theme_minimal() +
	  theme(
	    axis.text.x = element_text(angle = 45, hjust = 1),
	    axis.text.y = element_text(size = 10),
	    legend.title = element_text(size = 12),
	    legend.text = element_text(size = 10)
	  ) +
	  labs(
	    x = "Category",
	    y = "PRECAST_07",
	    title = "Heatmap of Proportions"
	  )
	dev.off()  
	
	spe_scz = readRDS(here("spe_scz.rds"))
	colData(spe_scz)$slide_id <- sapply(strsplit(colData(spe_scz)$sample_id, "_"), `[`, 1)
	df = as.data.frame(colData(spe_scz))
	
	
	# Create the heatmap by channel
	
	# Create a vector of unique PRECAST_07 values
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
	  
	  # Calculate the proportion for NeuN
	  proportion_neun <- sum(df$neun_pos == TRUE & df$PRECAST_07 == spd) / sum(df$neun_pos == TRUE)
	  prop_neun[[spd]] <- proportion_neun
	  
	  # Calculate the proportion for PNN
	  proportion_pnn <- sum(df$pnn_pos == TRUE & df$PRECAST_07 == spd) / sum(df$pnn_pos == TRUE)
	  prop_pnn[[spd]] <- proportion_pnn
	  
	  # Calculate the proportion for vascular position
	  proportion_vasc <- sum(df$vasc_pos == TRUE & df$PRECAST_07 == spd) / sum(df$vasc_pos == TRUE)
	  prop_vasc[[spd]] <- proportion_vasc
	}
	
	# Combine the results into a data frame
	prop_df1 <- data.frame(
	  PRECAST_07 = unique_spd_values,
	  Neuropil = unlist(prop_neuropil),
	  NeuN = unlist(prop_neun),
	  WFA = unlist(prop_pnn),
	  Claudin_5 = unlist(prop_vasc)
	)
	
		# Melt the dataframe for plotting
	prop_melted <- melt(prop_df, id.vars = "PRECAST_07") %>%
	  mutate(
	    PRECAST_07 = factor(PRECAST_07, levels = c("spd04","spd01","spd03", "spd05","spd02","spd06", "spd07"),
	                        labels = c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3","SpD07-L1" ))
	  )
	
	# Plot the heatmap
	plot1 <- ggplot(prop_melted[prop_melted$variable == "Neuropil",], aes(x = variable, y = PRECAST_07, fill = value)) +
	  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "black", limits = c(0,1))+  
	  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	  guides(fill = guide_legend(title = NULL)) + 
	  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
	plot2 <- ggplot(prop_melted[prop_melted$variable == "NeuN",], aes(x = variable, y = PRECAST_07, fill = value)) +
	  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "yellow", limits = c(0, 1)) + 
	  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	  guides(fill = guide_legend(title = NULL)) +
	  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
	plot3 <- ggplot(prop_melted[prop_melted$variable == "WFA",], aes(x = variable, y = PRECAST_07, fill = value)) +
	  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "magenta", limits = c(0, 1)) + 
	  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	  guides(fill = guide_legend(title = NULL)) +
	  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
	plot4 <- ggplot(prop_melted[prop_melted$variable == "Claudin_5",], aes(x = variable, y = PRECAST_07, fill = value)) +
	  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) + 
	  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	  guides(fill = guide_legend(title = NULL)) +
	  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3)
	
	plotC = grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
	ggsave("combined_plot.png", plot = plotC, width = 12, height = 10, dpi = 300)

	pdf(here("plots", "image_processing", "enrichment", "scz_combined.pdf"), width = 12, height = 10)
	ggplot(prop_melted, aes(x = variable, y = PRECAST_07, fill = value)) +
	  geom_tile(color = "white") + scale_fill_gradient(low = "white", high = "black", limits = c(0,1))+  
	  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	  guides(fill = guide_legend(title = NULL)) + 
	  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) 
	dev.off()
