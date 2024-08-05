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

gene <- rownames(spe)[match("PVALB", rowData(spe)$gene_name)]		  
spT <- spe[gene, ]
spT$gene = as.numeric(t(assays(spT)$counts))
df = as.data.frame(colData(spT))


df <- df %>%
  mutate(PNN = ifelse(PWFA > 0.05, "PNN+", "PNN-"))

# Calculate the proportion of PNN+ cells by PRECAST_07
prop_df <- df %>%
  group_by(PRECAST_07, dx) %>%
  summarize(proportion_PNN = mean(PNN == "PNN+"))

# Create the bar plot
ggplot(prop_df, aes(x = PRECAST_07, y = proportion_PNN, fill = dx)) +
  geom_bar(stat = "identity",  position = "dodge") +
  labs(title = "Proportion of PNN+ Cells by PRECAST_07",
       x = "PRECAST_07",
       y = "Proportion of PNN+ Cells") +
  theme_minimal()
  
  
  df <- df %>%
    mutate(PNN = ifelse(PWFA > 0.05, "PNN+", "PNN-"))

	prop_df <- df %>%
	  group_by(PRECAST_07, sample_id, dx) %>%
	  summarize(proportion_PNN = mean(PNN == "PNN+"))
	  
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
	