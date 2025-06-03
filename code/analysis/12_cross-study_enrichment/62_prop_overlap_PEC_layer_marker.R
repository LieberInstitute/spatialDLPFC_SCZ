# Load Library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
})


# Load marker genes df ----
## PEC Cell-type markers ----
cell_type_enrich_df <- readRDS(
  # NOTE: if run on JHPCE
  # file.path(
  #   # PEC Study folder
  #   "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/",
  #   "processed-data/rdata/spe/14_spatial_registration_PEC",
  #   "registration_stats_LIBD.rds"
  # )
  # NOTE: if run on Boyi's local machine
  here(
    "code/analysis/12_cross-study_enrichment",
    "registration_stats_LIBD.rds"
  )
)

cell_type_order <- c(
  # Excitatary neurons
  "L2.3.IT", "L4.IT", "L5.ET", "L5.IT",
  "L5.6.NP", "L6.CT", "L6.IT", "L6.IT.Car3", "L6b",
  # Inhibitory neurons
  "Chandelier", "Lamp5", "Lamp5.Lhx6", "Pax6", "Pvalb",
  "Sncg", "Sst", "Vip",
  # Non-neuronal cells
  "Astro", "Endo", "Immune", "Micro", "OPC", "Oligo",
  "PC", "SMC", "VLMC"
)

cell_tstats <- cell_type_enrich_df[, grep("[f|t]_stat_", colnames(cell_type_enrich_df))]
colnames(cell_tstats) <- gsub("[f|t]_stat_", "", colnames(cell_tstats))
cell_fdrs <- cell_type_enrich_df[, grep("fdr_", colnames(cell_type_enrich_df))]
colnames(cell_fdrs) <- gsub("fdr_", "", colnames(cell_fdrs))

### Create gene-by-layer DF of {0,1} ----
cell_marker_df <- cell_type_order |>
  set_names() |>
  map(
    ~ cell_tstats[, .x] > 0 & cell_fdrs[, .x] < 0.05,
  ) |>
  bind_cols()

colSums(cell_marker_df)


## Create heatmap of overlap proportion ----
cor(cell_marker_df)

cell_overlap_propo_matrix <- matrix(NA, ncol(cell_marker_df), ncol(cell_marker_df),
  dimnames = list(colnames(cell_marker_df), colnames(cell_marker_df))
)
for (i in seq_len(ncol(cell_marker_df))) {
  for (j in seq_len(ncol(cell_marker_df))) {
    cell_overlap_propo_matrix[i, j] <- sum(
      cell_marker_df[, i] & cell_marker_df[, j]
    ) / sum(
      cell_marker_df[, i] | cell_marker_df[, j]
    )
  }
}

# Visualize overlap proportion matrix ----
library(pheatmap)
pdf(
  here(
    "code/analysis/12_cross-study_enrichment",
    "prop_overlap_PEC_layer_marker.pdf"
  ),
  width = 10,
  height = 10
)

pheatmap(
  cor(cell_marker_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  # fontsize_row = 10,
  # fontsize_col = 10,
  main = "Correlation between PEC layer marker genes",
  color = colorRampPalette(c("red", "white", "blue"))(100),
  breaks = seq(-1, 1, length.out = 101)
)

pheatmap(
  cell_overlap_propo_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  # fontsize_row = 10,
  # fontsize_col = 10,
  main = "Proportion of overlap between PEC layer marker genes",
  color = colorRampPalette(c("white", "blue"))(100)
)
dev.off()
