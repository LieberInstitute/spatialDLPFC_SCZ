suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(spatialLIBD)
  library(ComplexHeatmap)
  library(circlize)
  library(sessioninfo)
})

# Load Data ----
## Load layer-adjusted DEGs ----
gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)

## Load PEC cell type ----
### Load Enrichment results with annotated cell types ----
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

#### Retrieve # of marker genes per cell type ----
tstats <- cell_type_enrich_df[, grep("[f|t]_stat_", colnames(cell_type_enrich_df))]
colnames(tstats) <- gsub("[f|t]_stat_", "", colnames(tstats))
fdrs <- cell_type_enrich_df[, grep("fdr_", colnames(cell_type_enrich_df))]
colnames(fdrs) <- gsub("fdr_", "", colnames(fdrs))

n_cell_type_marker <- colnames(tstats) |>
  set_names() |>
  map_dbl(
    # NOTE: the fdr threshold should be consistent with the enrichrment test parameters
    ~ sum(tstats[, .x] > 0 & fdrs[, .x] < 0.05)
  )


# Enrichment Analaysis ----
layer_adj_DEG_list <- list(
  # `Overall` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
  `Up` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
  `Down` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
)

overall_enrich_res <- spatialLIBD::gene_set_enrichment(
  gene_list = layer_adj_DEG_list,
  modeling_results = list(
    "enrichment" = cell_type_enrich_df
  ),
  model_type = "enrichment",
  fdr_cut = 0.05
)

# Enrichment Dot Plot ----
res <- overall_enrich_res

mat <- res |>
  mutate(
    OR = ifelse(OR < 1, 1, OR),
  ) |>
  select(
    ID, test, OR
  ) |>
  pivot_wider(names_from = test, values_from = OR, values_fill = 0) |>
  column_to_rownames("ID") |>
  as.matrix()

# Scale the size based on -log10(Pval)
size_mat <- res |>
  mutate(
    fdr_p = Pval |> p.adjust(method = "fdr"),
    sig_cat = case_when(
      Pval >= 0.05 ~ "nonsig",
      Pval < 0.05 & fdr_p >= 0.05 ~ "Nominal p < 0.05",
      fdr_p < 0.05 ~ "FDR < 0.05"
    ),
    size = case_when(
      sig_cat == "nonsig" ~ 0.5,
      sig_cat == "Nominal p < 0.05" ~ 1,
      sig_cat == "FDR < 0.05" ~ 1.5
    )
  ) |>
  select(
    ID, test, size
  ) |>
  pivot_wider(names_from = test, values_from = size, values_fill = 0) |>
  column_to_rownames("ID") |>
  as.matrix()

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

mat <- mat[, cell_type_order]
size_mat <- size_mat[, cell_type_order]

# browser()

celltype_ha <- HeatmapAnnotation(
  n_celltype = anno_barplot(
    n_cell_type_marker[cell_type_order],
    border = FALSE,
    axis = TRUE,
    gp = gpar(fill = "black"),
    height = unit(2, "cm"),
    axis_param = list(
      side = "left",
      at = seq(0, max(n_cell_type_marker), by = 1000),
      labels = seq(0, max(n_cell_type_marker), by = 1000)
    )
  ),
  gap = unit(5, "cm"),
  show_annotation_name = FALSE,
  name = "Huuki-Myers et al."
)

n_DEG <- res |>
  select(ID, SetSize) |>
  # group_by(ID) |>
  distinct(ID, SetSize) |>
  deframe()

DEG_ha <- rowAnnotation(
  n_DEG = anno_barplot(
    n_DEG,
    border = FALSE,
    axis = TRUE,
    gp = gpar(fill = c("Up" = "red", "Down" = "blue")),
    width = unit(1, "cm"),
    axis_param = list(
      direction = "reverse",
      at = seq(0, max(n_DEG), by = 50),
      labels = seq(0, max(n_DEG), by = 50),
      labels_rot = 0
    )
  ),
  gap = unit(5, "point"),
  show_annotation_name = FALSE

  # ,
  # annotation_name_side = "left",
  # annotation_name_gp = gpar(fontsize = 12)
)

# Change matrix orientation
# mat <- t(mat)
# size_mat <- t(size_mat)

# Define color function for Odds Ratio
col_fun <- colorRamp2(
  # TODO: change the color scale
  c(min(mat), median(mat), max(mat)),
  c("grey", "yellow", "blue")
)

# Create the dot plot using ComplexHeatmap
ht_list <- Heatmap(
  mat,
  name = "Odds Ratio",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  # Keep cell boundaries with black lines but hide the heatmap elements
  rect_gp = gpar(col = "lightgrey", lwd = 0.5, fill = NA),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x, y = y,
      r = unit(size_mat[i, j], "mm"),
      gp = gpar(fill = col_fun(mat[i, j]), col = NA)
    )
  },
  # Add cell type annotation - bar plot
  top_annotation = celltype_ha,
  left_annotation = DEG_ha,

  # Aesthetics
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8, rot = 45, just = "right"),
  heatmap_legend_param = list(
    title = "Odds Ratio",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 6) # ,
    # at = c(min(mat), max(mat))
  )
)


# draw(ht_list, heatmap_body = NULL)




# lgd_list <- list(
#   # dot size for p-value
#   Legend(
#     labels = c("Not sig.", "Nominal p < 0.05", "FDR < 0.05"),
#     title = "Significance", type = "points",
#     pch = 16,
#     legend_gp = gpar(fill = "black"),
#     size = unit(1:3, "mm"),
#   )
# )

# draw(ht_list,
#   annotation_legend_list = lgd_list,
#   column_title = title,
#   column_title_gp = grid::gpar(fontsize = 16)
# )

#   ht_list
# }


pdf(
  here(
    "plots/12_cross_study_enrichment",
    "layer_adj_DEG_vs_PEC_cell_type.pdf"
  ),
  height = 2, width = 5
)
ht_list |> draw(show_heatmap_legend = FALSE)
dev.off()

# Session Info ----
sessioninfo::session_info()
