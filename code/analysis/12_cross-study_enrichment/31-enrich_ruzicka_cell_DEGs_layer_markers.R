# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(spatialLIBD)
  library(sessioninfo)
})

# Load Data ----
## Load Ruzicka DEGs ----
# NOTE: clean up part in `code/analysis/12_cross-study_enrichment/clean_ruzicka_cell_DEG.R``
# Each element of the list is a data frame contianing cell-type specific DEGs
ruzicka_deg_list <- readRDS(here(
  "processed-data/rds/12_cross-study_enrichment",
  "ruzicka_cell_type_DEG.rds"
))

# Error prevention
stopifnot(
  length(ruzicka_deg_list) == 25
)

## Load PRECAST SpD09 layer-markers ----
raw_layer_df <- read_csv(
  here(
    "code/analysis/04_SpD_marker_genes",
    "raw_enrichment_spd_marker_gene.csv"
  )
)

## Load annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

# change the order of the labels
spd_order <- spd_anno_df$order

# Descriptive statistics ----
## all DEGs ----
ruzicka_deg_list_sig <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(ruzicka_sig_gene) |>
      pull(ensembl, name = NULL)
  )

ruzicka_deg_list_sig |>
  map_int(~ length(.x)) |>
  sum()
# [1] 6634

# Cell types with less than 25 DEGs
# Can't be used in spaitalLIBD enrichment analysis
which(ruzicka_deg_list_sig |>
  map_int(~ length(.x)) < 25)
#  Endo Pericytes


### Up-reg gene list ----
ruzicka_deg_list_up <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(
        ruzicka_sig_gene,
        Meta_logFC > 0
      ) |>
      pull(ensembl, name = NULL)
  )

# Number of up-reg DEGs per cell type
ruzicka_deg_list_up |>
  map_int(~ length(.x))
#           Ex-L2           Ex-L23            Ex-L3
#              135              186               67
#       Ex-L4_MYLK       Ex-L45_MET     Ex-L45_LRRK1
#               92               66              112
#     Ex-L5b_HTR2C           Ex-L56  Ex-L56_CC_NTNG2
#               56               54               29
#  Ex-L6_CC_SEMA3A    Ex-L6b_SEMA3D    Ex-L6b_SEMA3E
#               57               79               41
# In-Rosehip_CHST9 In-Rosehip_TRPC3        In-Reelin
#               23               57               50
#           In-VIP In-PV_Chandelier     In-PV_Basket
#               29               19               72
#           In-SST              Oli              OPC
#               19               42               44
#              Ast              Mic             Endo
#               83               70                1
#        Pericytes
#               10

### Down-reg gene list ----
ruzicka_deg_list_down <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(
        ruzicka_sig_gene,
        Meta_logFC < 0
      ) |>
      pull(ensembl, name = NULL)
  )
ruzicka_deg_list_down |>
  map_int(~ length(.x))
#            Ex-L2           Ex-L23            Ex-L3
#              251              362              375
#       Ex-L4_MYLK       Ex-L45_MET     Ex-L45_LRRK1
#              337              402              246
#     Ex-L5b_HTR2C           Ex-L56  Ex-L56_CC_NTNG2
#              226              348               54
#  Ex-L6_CC_SEMA3A    Ex-L6b_SEMA3D    Ex-L6b_SEMA3E
#              595              662              297
# In-Rosehip_CHST9 In-Rosehip_TRPC3        In-Reelin
#               86              241              145
#           In-VIP In-PV_Chandelier     In-PV_Basket
#              127               55              102
#           In-SST              Oli              OPC
#               75               21               34
#              Ast              Mic             Endo
#               65               26                6
#        Pericytes
#                3

# Run gene_set_enrichment_test ----

enrichment_dot_plot_ggplot <- function(
    res # , PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    ) {
  # browser()
  res |>
    ggplot(aes(x = test, y = ID, color = OR, size = -log10(Pval))) +
    geom_point() +
    scale_color_viridis_c(option = "C", direction = -1) + # Colorblind-friendly palette
    scale_size_continuous(range = c(1, 10)) +
    theme_minimal(base_size = 14) + # Adjust text size for publication
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    labs(
      color = "Odds Ratio",
      size = "-log10(p-value)"
    )
}


enrichment_dot_plot_heatmap <- function(
    res # , PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    ) {
  browser()
  library(ComplexHeatmap)
  library(circlize)

  # Prepare data for ComplexHeatmap
  mat <- res |>
    pivot_wider(names_from = test, values_from = OR, values_fill = 0) |>
    column_to_rownames("ID") |>
    as.matrix()

  # Scale the size based on -log10(Pval)
  size_mat <- res |>
    pivot_wider(names_from = test, values_from = -log10(Pval), values_fill = 0) |>
    column_to_rownames("ID") |>
    as.matrix()

  # Define color function for Odds Ratio
  col_fun <- colorRamp2(c(min(mat), max(mat)), c("white", "blue"))

  # Create the dot plot using ComplexHeatmap
  Heatmap(
    mat,
    name = "Odds Ratio",
    col = col_fun,
    rect_gp = gpar(type = "none"), # Remove default heatmap rectangles
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.circle(
        x = x, y = y,
        r = unit(size_mat[i, j] / max(size_mat), "mm"),
        gp = gpar(fill = col_fun(mat[i, j]), col = NA)
      )
    },
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12, rot = 45, just = "right"),
    heatmap_legend_param = list(
      title = "Odds Ratio",
      title_gp = gpar(fontsize = 14),
      labels_gp = gpar(fontsize = 12)
    )
  )
}

##  ALL GENS ----
all_gene_res <- ruzicka_deg_list_sig[
  # remove cell types with less than 25 DEGs
  -which(ruzicka_deg_list_sig |> map_int(~ length(.x)) < 25)
] |>
  spatialLIBD::gene_set_enrichment(
    gene_list = _,
    modeling_results = list("enrichment" = raw_layer_df),
    model_type = "enrichment",
    fdr_cut = 0.1
  ) |>
  # Annotate the spds
  inner_join(
    spd_anno_df,
    by = c("test" = "spd")
  ) |>
  mutate(test = anno_lab) |>
  select(-label, -anno_lab)

all_gene_res |> enrichment_dot_plot_ggplot()
all_gene_res |> enrichment_dot_plot_heatmap()

# |>
# gene_set_enrichment_plot(
#   # res,
#   PThresh = 12,
#   ORcut = 3,
#   enrichOnly = FALSE,
#   cex = 0.5 # control the size of the text
# ) + title(
#   "Ruzicka cell type-dx-DEGs enriched in PRECAST spd"
# )

# TODO: save the plot

## Up-reg gene only enrichment ----
# Remove cell ypes with less than 25 DEGs
ruzicka_deg_list_up <- ruzicka_deg_list_up[-which(ruzicka_deg_list_up |>
  map_int(~ length(.x)) < 25)]
ruzicka_deg_list_up |>
  spatialLIBD::gene_set_enrichment(
    modeling_results = list("enrichment" = raw_layer_df),
    model_type = "enrichment",
    fdr_cut = 0.1
  ) |>
  # Annotate the spds
  inner_join(
    spd_anno_df,
    by = c("test" = "spd")
  ) |>
  mutate(test = anno_lab) |>
  select(-label, -anno_lab) |>
  # make the plot
  gene_set_enrichment_plot(
    # res,
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    cex = 0.5 # control the size of the text
  ) + title(
    "Ruzicka cell type-dx-DEGs (up-reg only) enriched in PRECAST spd"
  )

# TODO: save the plot

## Down-reg gene only enrichment ----
ruzicka_deg_list_down <- ruzicka_deg_list_down[-which(ruzicka_deg_list_down |>
  map_int(~ length(.x)) < 25)]
ruzicka_deg_list_down |>
  spatialLIBD::gene_set_enrichment(
    modeling_results = list("enrichment" = raw_layer_df),
    model_type = "enrichment",
    fdr_cut = 0.1
  ) |>
  # Annotate the spds
  inner_join(
    spd_anno_df,
    by = c("test" = "spd")
  ) |>
  mutate(test = anno_lab) |>
  select(-label, -anno_lab) |>
  # make the plot
  gene_set_enrichment_plot(
    # res,
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    cex = 0.5 # control the size of the text
  ) + title(
    "Ruzicka cell type-dx-DEGs (down-reg only) enriched in PRECAST spd"
  )

# Session info ----
session_info()
