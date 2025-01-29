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

spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_w_spg_N63.rds"))

## Load SpD data ----
finalized_spd <- readRDS(here("processed-data/rds/spatial_cluster", "PRECAST","test_clus_label_df_semi_inform_k_2-16.rds"))

## Attach SpD label to spe ----
col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(spe$spg_PDAPI > 0.05,FALSE, TRUE)
spe$neun_pos <- ifelse(spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,TRUE, FALSE)
spe$vasc_pos <- ifelse(spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,TRUE, FALSE)

spe_ntc <- spe[, colData(spe)$dx == "ntc"]
spe_scz <- spe[, colData(spe)$dx == "scz"]

saveRDS(spe_ntc, here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
saveRDS(spe_scz, here("processed-data", "image_processing", "enrichment", "spe_scz.rds"))
####### ntc heatmap #########
spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
df = as.data.frame(colData(spe_ntc))

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
pdf(here("plots", "image_processing", "enrichment", "ntc.pdf"), width = 8, height = 6)
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


  ####### scz heatmap #########
  spe_scz = readRDS(here("processed-data", "image_processing", "enrichment", "spe_scz.rds"))
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
  pdf(here("plots", "image_processing", "enrichment", "scz.pdf"), width = 8, height = 6)
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