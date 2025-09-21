# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
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

ruzicka_cell_types <- names(ruzicka_deg_list)

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

length(ruzicka_deg_list_up)

# Number of up-reg DEGs per cell type
# ruzicka_deg_list_up |>
#   map_int(~ length(.x))
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

length(ruzicka_deg_list_down)

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
## Dot plot wrapper (complexHeatmap) ----
enrichment_dot_plot_heatmap <- function(
    res,
    #  PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    title) {
  up_flag <- stringr::str_detect(title, "Up")

  # Prepare data for ComplexHeatmap
  mat <- res |>
    mutate(
      OR = ifelse(OR < 1, 1, OR),
    ) |>
    select(
      ID, test, OR
    ) |>
    pivot_wider(names_from = test, values_from = OR, values_fill = 0) |>
    left_join(
      x = data.frame(ID = ruzicka_cell_types),
      by = "ID"
    ) |>
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
    left_join(
      x = data.frame(ID = ruzicka_cell_types),
      by = "ID"
    ) |>
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
  cell_type_order <- ruzicka_cell_types

  # TODO: manullay group cells to groups

  mat <- mat[cell_type_order, spd_order]
  size_mat <- size_mat[cell_type_order, spd_order]

  stopifnot(nrow(mat) == 25 && nrow(size_mat) == 25)



  # Change matrix orientation
  # mat <- t(mat)
  # size_mat <- t(size_mat)

  # Define color function for Odds Ratio
  col_fun <- colorRamp2(
    c(min(mat, na.rm = TRUE), median(mat, na.rm = TRUE), max(mat, na.rm = TRUE)),
    c("white", "yellow", "blue")
  )


  # Add row annotations to visualize number of DEGs per cell type

  deg_list <- list()

  if (up_flag) {
    deg_list <- ruzicka_deg_list_up
  } else {
    deg_list <- ruzicka_deg_list_down
  }

  if (length(deg_list) != nrow(mat)) {
    missing_names <- setdiff(rownames(mat), names(deg_list))
    if (length(missing_names) > 0) {
      for (nm in missing_names) {
        deg_list[[nm]] <- character(0)
      }
    }
    deg_list <- deg_list[rownames(mat)]
    stopifnot(length(deg_list) == nrow(mat))
  }

  row_anno <- rowAnnotation(
    Num_DEGs = anno_barplot(
      map_int(
        deg_list,
        ~ length(.x)
      )[rownames(mat)],
      gp = gpar(fill = ifelse(up_flag,
        "red", "blue"
      )),
      border = FALSE,
      axis = TRUE,
      axis_param = list(
        # at = c(0, 100, 200, 300),
        # labels = c(0, 100, 200, 300),
        gp = gpar(fontsize = 8),
        direction = ifelse(up_flag, "normal", "reverse"),
        labels_rot =  90 #ifelse(up_flag, -90, 90)
      ),
      width = unit(1, "cm"),
    ),
    show_annotation_name = FALSE # ,
    # annotation_name_gp = gpar(fontsize = 8),
    # annotation_name_side = "top"
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
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8, rot = 45, just = "right"),
    left_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Odds Ratio",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      at = c(1, 3, 6)
    )
  )

  return(ht_list)
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
  #   column_title_gp = gpar(fontsize = 8, fontface = "bold")
  # )
}

## Up-reg gene only enrichment ----
ruzicka_deg_list_up <- ruzicka_deg_list_up[
  # Remove cell ypes with less than 25 DEGs
  -which(ruzicka_deg_list_up |>
    map_int(~ length(.x)) < 25)
]
up_gene_res <- ruzicka_deg_list_up |>
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
  select(-label, -anno_lab)

pdf(
  here(
    "plots/12_cross_study_enrichment",
    "Ruzicka_DEGs_vs_PRECAST_layer_markers_up-reg.pdf"
  ),
  height = 4, width = 3
)
up_gene_res |>
  enrichment_dot_plot_heatmap(
    title = "Up-regulated Ruzicka DEGs"
  ) |>
  draw(show_heatmap_legend = FALSE)
dev.off()

## Down-reg gene only enrichment ----
ruzicka_deg_list_down <- ruzicka_deg_list_down[-which(ruzicka_deg_list_down |>
  map_int(~ length(.x)) < 25)]
down_gene_res <- ruzicka_deg_list_down |>
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
  select(-label, -anno_lab)

pdf(
  here(
    "plots/12_cross_study_enrichment",
    "Ruzicka_DEGs_vs_PRECAST_layer_markers_down-reg.pdf"
  ),
  height = 4, width = 3
)
down_gene_res |>
  enrichment_dot_plot_heatmap(
    title = "Down-regulated Ruzicka DEGs"
  ) |>
  draw(show_heatmap_legend = FALSE)
dev.off()

# Session info ----
session_info()


# Deprecated functions ----
##  ALL GENS ----
# all_gene_res <- ruzicka_deg_list_sig[
#   # remove cell types with less than 25 DEGs
#   -which(ruzicka_deg_list_sig |> map_int(~ length(.x)) < 25)
# ] |>
#   spatialLIBD::gene_set_enrichment(
#     gene_list = _,
#     modeling_results = list("enrichment" = raw_layer_df),
#     model_type = "enrichment",
#     fdr_cut = 0.1
#   ) |>
#   # Annotate the spds
#   inner_join(
#     spd_anno_df,
#     by = c("test" = "spd")
#   ) |>
#   mutate(test = anno_lab) |>
#   select(-label, -anno_lab)

# # all_gene_res |> enrichment_dot_plot_ggplot()
# pdf(
#   here(
#     "plots/12_cross_study_enrichment",
#     "Ruzicka_DEGs_vs_PRECAST_layer_markers_all.pdf"
#   ),
#   height = 4
# )
# all_gene_res |> enrichment_dot_plot_heatmap(
#   title = "Overall Ruzicka DEGs"
# )
# dev.off()

## (Deprecated) Dot plot wrapper (ggplot) ----
# enrichment_dot_plot_ggplot <- function(
#     res) {
#   res |>
#     ggplot(aes(x = test, y = ID, color = OR, size = -log10(Pval))) +
#     geom_point() +
#     # Colorblind-friendly palette
#     scale_color_viridis_c(option = "C", direction = -1) +
#     scale_size_continuous(range = c(1, 10)) +
#     theme_minimal(base_size = 14) + # Adjust text size for publication
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
#       axis.text.y = element_text(size = 12),
#       axis.title = element_blank(),
#       legend.title = element_text(size = 14),
#       legend.text = element_text(size = 12)
#     ) +
#     labs(
#       color = "Odds Ratio",
#       size = "-log10(p-value)"
#     )
# }
