# Load packages ----
# Load necessary packages and suppress startup messages
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(liana)
  library(ggplot2)
  library(sessioninfo)
})

# Load Data ----
## Load Spe Object that has spatial domain label ----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

PRECAST_df_final <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

col_data_df <- PRECAST_df_final |>
  right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

## SpD Annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))



# Run LIANA -----
## Prepare spe for LIANA format ----
spe$anno_spd <- spd_anno_df$anno_lab[match(spe$PRECAST_07, spd_anno_df$spd)]
colLabels(spe) <- spe$anno_spd


# change rownames(spe) to gene symbols
rownames(spe) <- rowData(spe)$gene_name

# Create dx-specific spe subset
.dx <- "ntc"
sub_spe <- spe[, spe$dx == .dx]
#sub_spe <- spe[, spe$sample_id == spe$sample_id[1]]

liana_test <- liana_wrap(
  sub_spe,
  seed = 20240925)

liana_res <- liana_test %>%
  liana_aggregate()

saveRDS(
  liana_res,
  here(
    "processed-data/layer_layer_comm",
    "test_liana_res_", .dx, ".pdf"
  )
)

# LIANA visualizations -----
pdf(
  here(
    "plots/layer_layer_comm",
  paste("test_liana_plots_" , .dx, ".pdf" )
  ),
  height = 12,
  width = 12
)


## Frequency Heatmap ----
liana_trunc <- liana_res %>%
   # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

heat_freq(liana_trunc)

## Circulized plot ----
p <- chord_freq(liana_trunc)
dev.off()

# Session Info ----
session_info()
