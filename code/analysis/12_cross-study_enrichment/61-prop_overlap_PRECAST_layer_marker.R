# Load Library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
})


# Load marker genes df ----
## PEC Cell-type markers ----
# cell_type_enrich_df <- readRDS(
#   # NOTE: if run on JHPCE
#   # file.path(
#   #   # PEC Study folder
#   #   "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/",
#   #   "processed-data/rdata/spe/14_spatial_registration_PEC",
#   #   "registration_stats_LIBD.rds"
#   # )
#   # NOTE: if run on Boyi's local machine
#   here(
#     "code/analysis/12_cross-study_enrichment",
#     "registration_stats_LIBD.rds"
#   )
# )


## PRECAST Spatial Domain Markers ----
PRECAST_layer_enrichment <- read_csv(
  here(
    "code/analysis/04_SpD_marker_genes",
    "raw_enrichment_spd_marker_gene.csv"
  )
)

## Load SpD annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = factor(
    paste0(gsub("spd", "SpD", spd), "-", label),
    levels = c(
      "SpD07-L1",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    )
  ))

### Number of markers per layer ----
layer_tstats <- PRECAST_layer_enrichment[, grep("[f|t]_stat_", colnames(PRECAST_layer_enrichment))]
colnames(layer_tstats) <- gsub("[f|t]_stat_", "", colnames(layer_tstats))
layer_fdrs <- PRECAST_layer_enrichment[, grep("fdr_", colnames(PRECAST_layer_enrichment))]
colnames(layer_fdrs) <- gsub("fdr_", "", colnames(layer_fdrs))


### Create gene-by-layer DF of {0,1} ----
layer_marker_df <- spd_anno_df$spd[order(spd_anno_df$anno_lab)] |>
  set_names(spd_anno_df$anno_lab[order(spd_anno_df$anno_lab)]) |>
  map(
    ~ layer_tstats[, .x] > 0 & layer_fdrs[, .x] < 0.05,
  ) |>
  bind_cols()

colSums(layer_marker_df)


# Calcualte the phi correaltion matrix ----
cor(layer_marker_df)

# Calcualte proportion of overlap -----
overlap_prop_matrix <- matrix(NA, ncol(layer_marker_df), ncol(layer_marker_df))
colnames(overlap_prop_matrix) <- rownames(overlap_prop_matrix) <- colnames(layer_marker_df)

for (i in seq_len(ncol(layer_marker_df))) {
  for (j in seq_len(ncol(layer_marker_df))) {
    n_overlap <- sum(layer_marker_df[[i]] & layer_marker_df[[j]])
    n_total_gene <- sum(
      layer_marker_df[, i] | layer_marker_df[, j]
    )
    overlap_prop_matrix[i, j] <- ifelse(n_total_gene == 0, NA, n_overlap / n_total_gene)
  }
}

# Visualize overlap_prop_matrix as a heatmap
library(pheatmap)



pdf(
  here(
    "code/analysis/12_cross-study_enrichment",
    "overlap_prop_matrix_layer_marker_genes.pdf"
  ),
  width = 10, height = 8
)

pheatmap(
  cor(layer_marker_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Correlation Between Layer marker genes",
  color = colorRampPalette(c("red", "white", "blue"))(100),
  breaks = seq(-1, 1, length.out = 101)
)


pheatmap(
  overlap_prop_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Proportion of Overlap Between Layer marker genes",
  color = colorRampPalette(c("white", "blue"))(100)
)
dev.off()
