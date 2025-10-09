# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(sessioninfo)
})

# Load data ----
## Load Ruzicka DEGs -----
ruzicka_deg_list <- readRDS(here(
  "processed-data/rds/12_cross-study_enrichment",
  "ruzicka_cell_type_DEG.rds"
))

ruzicka_cell_types <- names(ruzicka_deg_list)

ruzicka_deg_list_sig <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(ruzicka_sig_gene) |>
      pull(ensembl, name = NULL)
  )

ruzicka_deg_list_up <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(
        ruzicka_sig_gene,
        Meta_logFC > 0
      ) |>
      pull(ensembl, name = NULL)
  )

ruzicka_deg_list_down <- ruzicka_deg_list |>
  map(
    ~ .x |>
      filter(
        ruzicka_sig_gene,
        Meta_logFC < 0
      ) |>
      pull(ensembl, name = NULL)
  )

## Load Layer-restricted DEGs -----
layer_restricted_DEG_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
)

overall_deg_list <- layer_restricted_DEG_df |>
  filter(P.Value < 0.05) |>
  group_split(PRECAST_spd) |>
  set_names(
    nm = layer_restricted_DEG_df |>
      distinct(PRECAST_spd) |>
      pull(PRECAST_spd)
  )

up_reg_deg_list <- overall_deg_list |>
  map(~ .x |>
    filter(logFC > 0) |>
    pull(gene_id, name = NULL))

down_reg_deg_list <- overall_deg_list |>
  map(~ .x |>
    filter(logFC < 0) |>
    pull(gene_id, name = NULL))



## Create background gene sets ----
# Intersection of Ruzicka tested genes and layer-tested geness
ruzicka_deg_list |> map(~ .x |> nrow())

layer_restricted_DEG_df |>
  group_by(PRECAST_spd) |>
  summarize(
    n = n()
  )

## Create DEG gene sets -----



# Enrichment analysis ----

up_comp_df <- expand.grid(
  ruzicka_name = names(ruzicka_deg_list_up),
  layer_name = names(up_reg_deg_list),
  stringsAsFactors = FALSE
) |>
  pmap(
    .f = function(ruzicka_name, layer_name) {
      list(
        ruzicka_name = ruzicka_name,
        layer_name = layer_name,
        direction = "up",
        ruzicka_set = list(ruzicka_deg_list_up[[ruzicka_name]]),
        layer_set = list(up_reg_deg_list[[layer_name]]),
        background_set = list(intersect(
          ruzicka_deg_list[[ruzicka_name]] |> pull(ensembl, name = NULL),
          layer_restricted_DEG_df |>
            filter(PRECAST_spd == layer_name) |>
            pull(gene_id)
        ))
      )
    }
  ) |>
  bind_rows()

down_com_df <- expand.grid(
  ruzicka_name = names(ruzicka_deg_list_down),
  layer_name = names(down_reg_deg_list),
  stringsAsFactors = FALSE
) |>
  pmap_dfr(
    .f = function(ruzicka_name, layer_name) {
      list(
        ruzicka_name = ruzicka_name,
        layer_name = layer_name,
        direction = "down",
        ruzicka_set = list(ruzicka_deg_list_down[[ruzicka_name]]),
        layer_set = list(down_reg_deg_list[[layer_name]]),
        background_set = list(intersect(
          ruzicka_deg_list[[ruzicka_name]] |> pull(ensembl, name = NULL),
          layer_restricted_DEG_df |>
            filter(PRECAST_spd == layer_name) |>
            pull(gene_id)
        ))
      )
    }
  ) |>
  bind_rows()

## Fisher's exact test wrapper function ----
enrich_deg_aggreement <- function(
    set_1,
    set_2,
    n_background) {
  # browser()

  # create a list of contigency tables
  # set_1: vector of gene IDs (e.g., DEGs from Ruzicka)
  # set_2: vector of gene IDs (e.g., DEGs from layer-restricted)
  # n_background: total number of genes considered as background

  # Calculate overlap
  n_yes_yes <- length(intersect(set_1, set_2))
  n_yes_no <- length(set_1) - n_yes_yes
  n_no_yes <- length(set_2) - n_yes_yes
  n_no_no <- n_background - (n_yes_yes + n_yes_no + n_no_yes)

  # Build 2x2 contingency table
  cont_tab <- matrix(
    c(n_yes_yes, n_yes_no, n_no_yes, n_no_no),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      "set_1" = c("yes", "no"),
      "set_2" = c("yes", "no")
    )
  )

  fet_res <- fisher.test(cont_tab, alternative = "greater")

  list(
    OR = fet_res$estimate,
    Pval = fet_res$p.value,
    n_overlap = n_yes_yes,
    n_set_1 = length(set_1),
    n_set_2 = length(set_2),
    n_background = n_background
  )
}



# Run wrapper functions ----
up_enrich_long <- up_comp_df |>
  pmap(.f = function(ruzicka_name, layer_name,
                     ruzicka_set, layer_set, background_set,
                     ...) {
    ret_raw <- enrich_deg_aggreement(
      set_1 = ruzicka_set,
      set_2 = layer_set,
      n_background = length(background_set)
    )

    ret <- c(
      ruzicka_name = ruzicka_name,
      layer_name = layer_name,
      direction = "up",
      ret_raw
    )

    return(ret)
  }) |>
  bind_rows()

down_enrich_long <- down_com_df |>
  pmap(.f = function(ruzicka_name, layer_name,
                     ruzicka_set, layer_set, background_set,
                     ...) {
    ret_raw <- enrich_deg_aggreement(
      set_1 = ruzicka_set,
      set_2 = layer_set,
      n_background = length(background_set)
    )

    ret <- c(
      ruzicka_name = ruzicka_name,
      layer_name = layer_name,
      direction = "down",
      ret_raw
    )

    return(ret)
  }) |>
  bind_rows()

# Plot results ----

## Wrapper function -----
enrichment_dot_plot_heatmap <- function(
    res,
    #  PThresh = 12, ORcut = 3, enrichOnly = FALSE, cex = 0.5
    title) {
  # browser()
  up_flag <- stringr::str_detect(title, "Up")

  # Prepare data for ComplexHeatmap
  mat <- res |>
    mutate(
      OR = ifelse(OR < 1, 1, OR),
    ) |>
    select(
      ruzicka_name,
      test = layer_name, OR
    ) |>
    pivot_wider(names_from = test, values_from = OR, values_fill = 0) |>
    left_join(
      x = data.frame(ID = ruzicka_cell_types),
      by = c("ID" = "ruzicka_name")
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
      ID = ruzicka_name, test = layer_name, size
    ) |>
    pivot_wider(
      names_from = test, values_from = size, values_fill = 0
    ) |>
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
        labels_rot = 90 # ifelse(up_flag, -90, 90)
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


## Up-regulated enrichment ----
plot_up_enrich <- up_enrich_long |>
  mutate(
    # ruzicka_name = factor(ruzicka_name, levels = names(ruzicka_deg_list_up)),
    layer_name = factor(layer_name,
      levels = c(
        "SpD07-L1",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      ) |> rev(),
      labels = c(
        "SpD07-L1/M",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      ) |> rev()
    )
  ) |>
  enrichment_dot_plot_heatmap(
    title = "Ruzicka DEGs vs Layer-restricted DEGs (Up-regulated)"
  )

pdf(
  here(
    "plots/12_cross_study_enrichment",
    "ruzicka_DEG_vs_layer_DEGs_up.pdf"
  ),
  height = 4, width = 3
)
plot_up_enrich |>
  draw(show_heatmap_legend = FALSE)
dev.off()

# ggplot(aes(
#   x = ruzicka_name,
#   y = layer_name,
#   color = OR,
#   size = -log10(Pval)
# )) +
# geom_point() +
# scale_color_gradient2(
#   low = "blue", mid = "white", high = "red",
#   midpoint = 1, limits = c(0.5, max(up_enrich_long$OR, na.rm = TRUE)),
#   na.value = "lightgray"
# ) +
# labs(
#   title = "Ruzicka DEGs vs Layer-restricted DEGs (Up-regulated)",
#   x = "Ruzicka DEGs",
#   y = "Layer-restricted DEGs",
#   color = "Odds Ratio",
#   size = "-log10(P-value)"
# ) +
# theme_minimal() +
# theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_down_enrich <- down_enrich_long |>
  mutate(
    ruzicka_name = factor(ruzicka_name, levels = names(ruzicka_deg_list_down)),
    layer_name = factor(layer_name,
      levels = c(
        "SpD07-L1",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      ) |> rev(),
      labels = c(
        "SpD07-L1/M",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      ) |> rev()
    )
  ) |>
  enrichment_dot_plot_heatmap(
    title = "Ruzicka DEGs vs Layer-restricted DEGs (Down-regulated)"
  )

pdf(
  here(
    "plots/12_cross_study_enrichment",
    "ruzicka_DEG_vs_layer_DEGs_down.pdf"
  ),
  height = 4, width = 3
)
plot_down_enrich |>
  draw(show_heatmap_legend = FALSE)
dev.off()
## Save plots ----
# ggsave(
#   filename = here(
#     "plots/12_cross_study_enrichment",
#     "ruzicka_DEG_vs_layer_DEGs_up.pdf"
#   ),
#   plot = plot_up_enrich,
#   width = 8, height = 3
# )
# ggsave(
#   filename = here(
#     "plots/12_cross_study_enrichment",
#     "ruzicka_DEG_vs_layer_DEGs_down.pdf"
#   ),
#   plot = plot_down_enrich,
#   width = 8, height = 3
# )


# Session Info ----
session_info()
