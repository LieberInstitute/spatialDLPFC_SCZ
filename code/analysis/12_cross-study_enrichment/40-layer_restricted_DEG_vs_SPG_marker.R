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
## Load layer-restricted DEGs ----
layer_specifc_DEG_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
) |>
  # keep the nominal p-value < 0.05 DEGs
  filter(P.Value < 0.05)

overall_deg_list <- layer_specifc_DEG_df |>
  group_split(PRECAST_spd) |>
  set_names(
    nm = layer_specifc_DEG_df |>
      distinct(PRECAST_spd) |>
      pull(PRECAST_spd)
  )

up_reg_deg_list <- overall_deg_list |>
  map(~ .x |> filter(logFC > 0))

down_reg_deg_list <- overall_deg_list |>
  map(~ .x |> filter(logFC < 0))

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

#### Number of markers per SPG ----
n_spg_marker <- list(
  Neuropil = neuropil_df |> filter(fdr_TRUE < 0.05 & t_stat_TRUE > 0) |> nrow(),
  Neun = neun_df |> filter(fdr_TRUE < 0.05 & t_stat_TRUE > 0) |> nrow(),
  Vasculature = vasc_df |> filter(fdr_TRUE < 0.05 & t_stat_TRUE > 0) |> nrow(),
  PNN = pnn_df |> filter(fdr_TRUE < 0.05 & t_stat_TRUE > 0) |> nrow()
) |> unlist()



# Exact test for enrichment ----
## Overall DEGs ----
### Neuropil ----
neuropil_enrich_overall <- spatialLIBD::gene_set_enrichment(
  gene_list = overall_deg_list |> map(~ .x |> pull(gene_id)),
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

### Neun ----
neun_enrich_overall <- spatialLIBD::gene_set_enrichment(
  gene_list = overall_deg_list |> map(~ .x |> pull(gene_id)),
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

### Vasculature ----
vasc_enrich_overall <- spatialLIBD::gene_set_enrichment(
  gene_list = overall_deg_list |> map(~ .x |> pull(gene_id)),
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

### PNN ----
pnn_enrich_overall <- spatialLIBD::gene_set_enrichment(
  gene_list = overall_deg_list |> map(~ .x |> pull(gene_id)),
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

### Combine results ----
enrich_df_overall <- dplyr::bind_rows(
  neuropil_enrich_overall,
  neun_enrich_overall,
  vasc_enrich_overall,
  pnn_enrich_overall
) |> mutate(ID = gsub("SpD07-L1", "SpD07-L1/M", ID))

## Up-regulated DEGs ----
### Neuropil ----
neuropil_enrich_up <- spatialLIBD::gene_set_enrichment(
  gene_list = up_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### Neun ----
neun_enrich_up <- spatialLIBD::gene_set_enrichment(
  gene_list = up_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### Vasculature ----
vasc_enrich_up <- spatialLIBD::gene_set_enrichment(
  gene_list = up_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### PNN ----
pnn_enrich_up <- spatialLIBD::gene_set_enrichment(
  gene_list = up_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### Combine results ----
enrich_df_up <- dplyr::bind_rows(
  neuropil_enrich_up,
  neun_enrich_up,
  vasc_enrich_up,
  pnn_enrich_up
) |> mutate(ID = gsub("SpD07-L1", "SpD07-L1/M", ID))

## Down-regulated DEGs ----
### Neuropil ----
neuropil_enrich_down <- spatialLIBD::gene_set_enrichment(
  gene_list = down_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### Neun ----
neun_enrich_down <- spatialLIBD::gene_set_enrichment(
  gene_list = down_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### Vasculature ----
vasc_enrich_down <- spatialLIBD::gene_set_enrichment(
  gene_list = down_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### PNN ----
pnn_enrich_down <- spatialLIBD::gene_set_enrichment(
  gene_list = down_reg_deg_list |> map(~ .x |> pull(gene_id)),
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

### Combine results ----
enrich_df_down <- dplyr::bind_rows(
  neuropil_enrich_down,
  neun_enrich_down,
  vasc_enrich_down,
  pnn_enrich_down
) |> mutate(ID = gsub("SpD07-L1", "SpD07-L1/M", ID))

# Make dot plots ----
## Wrapper Function ----
enrichment_dot_plot_heatmap <- function(
    res,
    # PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    title) {
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

  if_order <- res$test |> unique()
  # browser()
  mat <- mat[spd_order, if_order]
  size_mat <- size_mat[spd_order, if_order]

  # browser()

  # spg_ha <- NULL
  # if (stringr::str_detect(title, "Down")) {
    spg_ha <- HeatmapAnnotation(
      spg = anno_barplot(
        n_spg_marker[if_order],
        axis = TRUE,
        axis_name = "# of SPG",
        border = FALSE,
        gp = gpar(fill = "black"),
        axis_param = list(
          at = seq(0, max(n_spg_marker[if_order]), by = 2000),
          labels = seq(0, max(n_spg_marker[if_order]), by = 2000),
          # gp = gpar(fontsize = 8),
          direction = "reverse",
          side = "right"
        ),
        Height = unit(2, "cm"),
      ),
      show_annotation_name = FALSE
    )
  # }

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
    # Add the SPG annotation
    bottom_annotation = spg_ha,
    # Aethetics
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8, rot = 45, just = "right"),
    heatmap_legend_param = list(
      title = "Odds Ratio",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8) # ,
      # at = c(min(mat), max(mat))
    )
  )

  return(ht_list)
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

  # draw(ht_list,
  #   # annotation_legend_list = lgd_list,
  #   column_title = title,
  #   column_title_gp = grid::gpar(fontsize = 8)
  # )
}

replace_Layer1_label <- function(x) {
  # browser()
  x <- gsub("SpD07-L1", "SpD07-L1/M", x)
  return(x)
}


## Up-regulated DEGs ----
pdf(
  here(
    "plots/12_cross_study_enrichment",
    "layer_restricted_DEG_vs_SPG_marker_upregulated.pdf"
  ),
  height = 2.3, width = 1.5
)
enrichment_dot_plot_heatmap(
  res = enrich_df_up,
  title = "Up-regulated DEGs"
) |> draw(show_heatmap_legend = FALSE)
dev.off()

## Down-regulated DEGs ----
pdf(
  here(
    "plots/12_cross_study_enrichment",
    "layer_restricted_DEG_vs_SPG_marker_downregulated.pdf"
  ),
  height = 2.3, width = 1.5
)
enrichment_dot_plot_heatmap(
  res = enrich_df_down,
  title = "Down-regulated DEGs"
) |> draw(show_heatmap_legend = FALSE)
dev.off()

# Session info ----
session_info()

# Deprecated code ----
## Overall DEGs ----
# pdf(
#   here(
#     "plots/12_cross_study_enrichment",
#     "layer_restricted_DEG_vs_SPG_marker_overall.pdf"
#   ),
#   height = 2.7, width = 2
# )
# enrichment_dot_plot_heatmap(
#   res = enrich_df_overall,
#   title = "Overall DEGs"
# )
# dev.off()
