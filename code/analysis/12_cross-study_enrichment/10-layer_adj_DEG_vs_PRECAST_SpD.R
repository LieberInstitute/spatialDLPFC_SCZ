# Load library ----
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

## Load layer-markers ----
PRECAST_layer_enrichment <- read_csv(
  here(
    "code/analysis/04_SpD_marker_genes",
    "raw_enrichment_spd_marker_gene.csv"
  )
)

### Number of markers per layer ----
tstats <- PRECAST_layer_enrichment[, grep("[f|t]_stat_", colnames(PRECAST_layer_enrichment))]
colnames(tstats) <- gsub("[f|t]_stat_", "", colnames(tstats))
fdrs <- PRECAST_layer_enrichment[, grep("fdr_", colnames(PRECAST_layer_enrichment))]
colnames(fdrs) <- gsub("fdr_", "", colnames(fdrs))

n_layer_marker <- colnames(tstats) |>
  set_names() |>
  map_dbl(
    # NOTE: the fdr threshold should be consistent with the enrichrment test parameters
    ~ sum(tstats[, .x] > 0 & fdrs[, .x] < 0.05)
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
      "SpD07-L1/M",
      "SpD06-L2/3",
      "SpD02-L3/4",
      "SpD05-L5",
      "SpD03-L6",
      "SpD01-WMtz",
      "SpD04-WM"
    )
  ))

# Enrichment test ----
layer_adj_DEG_list <- list(
  `Overall` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
  `Up` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
  `Down` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
)

res <- gene_set_enrichment(
  gene_list = layer_adj_DEG_list,
  modeling_results = list("enrichment" = PRECAST_layer_enrichment),
  model_type = "enrichment",
  fdr_cut = 0.05
)

# Make visualizaiton of the results ----
## Annotaed SpDs to layers ----
annotated_res <- res |>
  inner_join(
    spd_anno_df,
    by = c("test" = "spd")
  ) |>
  mutate(test = anno_lab) |>
  select(-label, -anno_lab)

n_layer_marker_annotated <- n_layer_marker
names(n_layer_marker_annotated) <- spd_anno_df$anno_lab[match(names(n_layer_marker), spd_anno_df$spd)]



## Visaulize via dot plot -----
# Wrapper function for the dot plot
enrichment_dot_plot_heatmap <- function(
    res # , PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    ) {
  # Prepare data for ComplexHeatmap
  # browser()
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

  # Reorder the matrix based on the order of spd_anno_df
  spd_order <- c(
    "SpD07-L1/M",
    "SpD06-L2/3",
    "SpD02-L3/4",
    "SpD05-L5",
    "SpD03-L6",
    "SpD01-WMtz",
    "SpD04-WM"
  )

  cell_type_order <- res$ID |> unique()

  mat <- mat[cell_type_order, spd_order]
  size_mat <- size_mat[cell_type_order, spd_order]

  # Change matrix orientation
  # Change to layer-by-DEGs
  mat <- t(mat)
  size_mat <- t(size_mat)

  # browser()

  n_degs <- res |>
    select(ID, SetSize) |>
    distinct(ID, SetSize) |>
    deframe()

  deg_ha <- HeatmapAnnotation(
    # Create a named list for the annotation
    DEGs = anno_barplot(
      n_degs,
      border = TRUE,
      axis = TRUE,
      gp = gpar(fill = "black"),
      height = unit(1, "cm"),
      axis_param = list(
        # direction = "reverse",
        at = seq(0, max(n_degs), by = 50),
        labels = seq(0, max(n_degs), by = 50)
      )
    ),
    annotation_name_gp = gpar(fontsize = 8) # set annotation title font size
  )

  layer_marker_ha <- rowAnnotation(
    # Create a named list for the annotation
    n_layer_marker = anno_barplot(
      n_layer_marker_annotated[spd_order],
      border = TRUE,
      axis = TRUE,
      gp = gpar(fill = "black"),
      height = unit(1, "cm"),
      axis_param = list(
        direction = "reverse",
        at = seq(0, max(n_layer_marker_annotated), by = 2000),
        labels = seq(0, max(n_layer_marker_annotated), by = 2000)
      )
    ),
    annotation_name_gp = gpar(fontsize = 8) # set annotation title font size
  )



  # Define color function for Odds Ratio
  col_fun <- colorRamp2(
    # TODO: change the color scale
    c(min(mat), median(mat), max(mat)),
    c("grey", "yellow", "blue")
  )

  # browser()

  # Create the dot plot using ComplexHeatmap
  ht_list <- Heatmap(
    mat,
    name = "Odds Ratio",
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    # Keep cell boundaries with black lines but hide the heatmap elements
    rect_gp = gpar(col = "black", lwd = 0.5, fill = NA),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.circle(
        x = x, y = y,
        r = unit(size_mat[i, j], "mm"),
        gp = gpar(fill = col_fun(mat[i, j]), col = NA)
      )
    },
    # Annotation Bar plot
    top_annotation = deg_ha,
    left_annotation = layer_marker_ha,

    # Aesthetics
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8, rot = 45, just = "right"),
    heatmap_legend_param = list(
      title = "Odds Ratio",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 6),
      at = c(1, 3, 6)
    )
  )

  # browser()
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

  # draw(ht_list, annotation_legend_list = lgd_list)

  ht_list
}


# Create complexHeatmap
ret_htmap <- enrichment_dot_plot_heatmap(
  res = annotated_res
)


# Save plot
pdf(
  file = here(
    "plots/12_cross_study_enrichment",
    "dotplot_layer_adj_DEG_vs_PRECAST_SpD.pdf"
  ),
  width = 3, height = 2.5
)
ret_htmap
dev.off()


# Session Info ----
sessioninfo::session_info()

# Deprecateted code
## Visaulize via heatmap (Deprecated) -----
# pdf(
#   file = here(
#     "plots/12_cross_study_enrichment",
#     "heatmap_layer_adj_DEG_vs_PRECAST_SpD.pdf"
#   ),
#   width = 3, height = 4
# )
# # plot
# annotated_res |>
#   gene_set_enrichment_plot(
#     # res,
#     PThresh = 12,
#     ORcut = 3,
#     enrichOnly = FALSE,
#     cex = 1.5 # control the size of the text
#   ) + title("Layer-adjusted SCZ-DEGs")
# dev.off()
