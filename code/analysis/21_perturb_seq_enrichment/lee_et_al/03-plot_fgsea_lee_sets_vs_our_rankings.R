# Plot fgsea results: Lee et al. perturbation gene sets against our rankings ----
# - Check against their Fig 6C (Allcells Matr3 and Scyl1)
# - NES heatmaps per family of our rankings: rows are their gene sets
#   (cell-type grouping x perturbation), columns are our rankings; down and up
#   sets side by side. Asterisks mark padj_family < 0.05. Crossed-out cells
#   were not tested (fewer than minSize genes after universe matching).
# - Bars: row bars are their DEG set sizes, column bars are our ranked genes,
#   each split by one-to-one ortholog status (see utils_heatmap_bars.R).

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(sessioninfo)
})
source(here("code/analysis/21_perturb_seq_enrichment/lee_et_al/utils_heatmap_bars.R"))


# Config ----
fdr_cut <- 0.05

fld_in <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")
fld_plot <- here("plots/21_perturb_seq_enrichment/lee_et_al")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)

celltype_order <- c(
  "Allcells", "Neuron",
  "Glutamatergic", "GABAergic", "Astro", "Oligo", "Microglia", "Vascular"
)

# Diverging NES scale: blue (down with SCZ / PRS) - gray - red (up)
nes_col <- colorRamp2(c(-3, 0, 3), c("#2a78d6", "#f0efec", "#e34948"))
untested_col <- "white"


# Load data ----
res_df <- readRDS(file.path(fld_in, "fgsea_lee_sets_vs_our_rankings.rds")) |>
  mutate(sig = if_else(!is.na(padj_family) & padj_family < fdr_cut, "*", ""))
ranking_meta <- readRDS(file.path(fld_in, "our_rankings.rds"))$meta
bd <- load_bar_data(fld_in)

pert_order <- sort(unique(res_df$pert))
row_df <- expand_grid(celltype = celltype_order, pert = pert_order) |>
  mutate(row_id = paste(celltype, pert, sep = "|"))


# Helper ----
# NES, tested and label matrices for one direction:
# rows = celltype x perturbation, columns = rankings
build_mats <- function(direction, rankings) {
  df <- res_df |>
    filter(dir == direction, ranking %in% rankings) |>
    mutate(row_id = paste(celltype, pert, sep = "|"))

  to_mat <- function(value_col, fill) {
    df |>
      select(row_id, ranking, value = all_of(value_col)) |>
      complete(row_id = row_df$row_id, ranking = rankings, fill = list(value = fill)) |>
      pivot_wider(names_from = ranking, values_from = value) |>
      column_to_rownames("row_id") |>
      as.matrix() |>
      (\(m) m[row_df$row_id, rankings, drop = FALSE])()
  }

  list(
    nes = to_mat("NES", NA_real_),
    tested = to_mat("tested", FALSE),
    sig = to_mat("sig", "")
  )
}

# Cell function: cross out untested cells, star significant ones
sig_cell_fun <- function(tested, sig) {
  function(j, i, x, y, w, h, fill) {
    if (!tested[i, j]) {
      grid.lines(
        x = unit.c(x - w / 2, x + w / 2),
        y = unit.c(y - h / 2, y + h / 2),
        gp = gpar(col = "#c3c2b7", lwd = 0.5)
      )
    } else if (sig[i, j] != "") {
      grid.text(sig[i, j], x, y, vjust = 0.75, gp = gpar(fontsize = 8))
    }
  }
}

make_ht <- function(direction, rankings, show_legend) {
  mats <- build_mats(direction, rankings)
  dir_label <- str_to_title(direction)

  Heatmap(
    mats$nes,
    name = paste0("NES_", direction),
    col = nes_col,
    na_col = untested_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(row_df$celltype, levels = celltype_order),
    row_labels = row_df$pert,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8),
    row_gap = unit(1.5, "mm"),
    column_title = paste(dir_label, "DEGs"),
    column_title_gp = gpar(fontsize = 9),
    column_labels = display_label(rankings),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 90,
    rect_gp = gpar(col = "#fcfcfb", lwd = 1),
    left_annotation = bar_anno(
      lee_set_counts(bd, row_df$celltype, row_df$pert, direction),
      "row",
      paste(dir_label, "DEGs")
    ),
    top_annotation = bar_anno(
      our_ranking_counts(bd, rankings),
      "column",
      "Ranked genes",
      size = unit(1.2, "cm")
    ),
    show_heatmap_legend = show_legend,
    heatmap_legend_param = list(title = "NES", at = c(-3, -1.5, 0, 1.5, 3)),
    cell_fun = sig_cell_fun(mats$tested, mats$sig)
  )
}

annotation_legends <- function() {
  list(
    Legend(
      labels = c(
        paste0("padj_family < ", fdr_cut, " (*)"),
        "Not tested (< minSize)"
      ),
      title = "",
      type = "points",
      pch = c("*", "/"),
      legend_gp = gpar(col = c("black", "#c3c2b7"))
    ),
    bar_legend()
  )
}

plot_family <- function(rankings, title, file, width) {
  ht <- make_ht("down", rankings, TRUE) + make_ht("up", rankings, FALSE)

  pdf(file.path(fld_plot, file), width = width, height = 17)
  draw(
    ht,
    column_title = title,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    ht_gap = unit(14, "mm"),
    annotation_legend_list = annotation_legends()
  )
  dev.off()
}


# Check against Lee et al. Fig 6C ----
fig6c_rankings <- c(
  "PRS",
  ranking_meta$ranking[ranking_meta$family == "per_domain"]
)
fig6c_rows <- expand_grid(pert = c("Matr3", "Scyl1"), dir = c("down", "up")) |>
  mutate(row_id = paste(pert, dir))

fig6c_df <- res_df |>
  filter(
    celltype == "Allcells",
    pert %in% c("Matr3", "Scyl1"),
    ranking %in% fig6c_rankings
  ) |>
  mutate(row_id = paste(pert, dir))

fig6c_mat <- function(value_col) {
  fig6c_df |>
    select(row_id, ranking, value = all_of(value_col)) |>
    pivot_wider(names_from = ranking, values_from = value) |>
    column_to_rownames("row_id") |>
    as.matrix() |>
    (\(m) m[fig6c_rows$row_id, fig6c_rankings, drop = FALSE])()
}
fig6c_nes <- fig6c_mat("NES")
fig6c_sig <- fig6c_mat("sig")

pdf(file.path(fld_plot, "fig6C_check_Allcells_Matr3_Scyl1.pdf"), width = 7.5, height = 3.6)
draw(
  Heatmap(
    fig6c_nes,
    name = "NES",
    col = nes_col,
    na_col = untested_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = factor(
      if_else(fig6c_rankings == "PRS", "PRS", "Layer-restricted SCZ"),
      levels = c("PRS", "Layer-restricted SCZ")
    ),
    column_gap = unit(2, "mm"),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontsize = 9),
    rect_gp = gpar(col = "#fcfcfb", lwd = 1),
    left_annotation = bar_anno(
      lee_set_counts(bd, "Allcells", fig6c_rows$pert, fig6c_rows$dir),
      "row", "DEGs"
    ),
    top_annotation = bar_anno(
      our_ranking_counts(bd, fig6c_rankings),
      "column", "Ranked genes", size = unit(1.2, "cm")
    ),
    heatmap_legend_param = list(title = "NES", at = c(-3, -1.5, 0, 1.5, 3)),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(
        sprintf("%.1f%s", fig6c_nes[i, j], fig6c_sig[i, j]),
        x, y,
        gp = gpar(fontsize = 7)
      )
    }
  ),
  column_title = "Allcells Matr3 / Scyl1 DEGs (cf. Lee et al. Fig 6C)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_list = list(bar_legend())
)
dev.off()


# Family heatmaps ----
plot_family(
  rankings = c("Domain-adjusted", "PRS"),
  title = "Layer-adjusted SCZ and PRS rankings",
  file = "heatmap_NES_domain_adjusted_PRS.pdf",
  width = 7
)

plot_family(
  rankings = ranking_meta$ranking[ranking_meta$family == "per_domain"],
  title = "Layer-restricted SCZ rankings",
  file = "heatmap_NES_per_domain.pdf",
  width = 10
)

plot_family(
  rankings = ranking_meta$ranking[ranking_meta$family == "spg"],
  title = "SPG-specific SCZ rankings",
  file = "heatmap_NES_spg.pdf",
  width = 8
)


# Session info ----
sessioninfo::session_info()
