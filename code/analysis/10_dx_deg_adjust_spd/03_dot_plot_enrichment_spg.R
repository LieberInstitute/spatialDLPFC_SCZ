# Load Library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  library(circlize)
  library(ComplexHeatmap)
  library(here)
  library(sessioninfo)
})

# Load Data ----
## Load DEGs ----
gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)

## Load SPG enrichment res ----
### Neuropil ----
neuropil_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "neuropil_dx_res.csv"
  )
)

### Neun ----
neun_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "neun_dx_res.csv"
  )
)

### Vasculature ----
vasc_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "vasc_dx_res.csv"
  )
)

### pnn ----
pnn_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "pnn_dx_res.csv"
  )
)

# Exact test for enrichment ----
## Neuropil ----
neuropil_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = list("enrichment" = neuropil_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  # Finding the marker genes
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "Neuropil"
  )

## Neun ----
neun_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = list("enrichment" = neun_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "Neun"
  )

## Vasculature ----
vasc_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = list("enrichment" = vasc_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "Vasculature"
  )

## PNN ----
pnn_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = list("enrichment" = pnn_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "PNN"
  )


# Combine results ----
enrich_df <- dplyr::bind_rows(
  neuropil_enrich,
  neun_enrich,
  vasc_enrich,
  pnn_enrich
)


## Make the dot plot ----
# TODO: copy the code here
enrichment_dot_plot_heatmap <- function(
    res # , PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    ) {
  # Prepare data for ComplexHeatmap
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
  # spd_order <- c(
  #   "SpD07-L1/M",
  #   "SpD06-L2/3",
  #   "SpD02-L3/4",
  #   "SpD05-L5",
  #   "SpD03-L6",
  #   "SpD01-WMtz",
  #   "SpD04-WM"
  # )

  # cell_type_order <- res$ID |> unique()

  # mat <- mat[cell_type_order, spd_order]
  # size_mat <- size_mat[cell_type_order, spd_order]

  # Change matrix orientation
  # mat <- t(mat)
  # size_mat <- t(size_mat)

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
        r = unit(1.5 * size_mat[i, j], "mm"),
        gp = gpar(fill = col_fun(mat[i, j]), col = NA)
      )
    },
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12, rot = 45, just = "right") # ,
    # heatmap_legend_param = list(
    #   title = "Odds Ratio",
    #   title_gp = gpar(fontsize = 14),
    #   labels_gp = gpar(fontsize = 12),
    #   at = c(1, 3, 6)
    # )
  )

  # browser()
  lgd_list <- list(
    # dot size for p-value
    Legend(
      labels = c("Not sig.", "Nominal p < 0.05", "FDR < 0.05"),
      title = "Significance", type = "points",
      pch = 16,
      legend_gp = gpar(fill = "black"),
      size = unit(1:3, "mm"),
    )
  )

  draw(ht_list, annotation_legend_list = lgd_list)
}

pdf(
  here(
    "plots/10_dx_deg_adjust_spd",
    "dotplot_layer_adj_DEG_enrich_PRECAST_markers.pdf"
  ),
  width = 4.5,
  height = 2
)
enrich_df |> enrichment_dot_plot_heatmap()
dev.off()

# Session info ----
session_info()