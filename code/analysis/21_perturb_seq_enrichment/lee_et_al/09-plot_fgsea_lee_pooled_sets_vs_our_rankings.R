# Plot fgsea results: Lee et al. pooled (up + down) DEG sets vs our rankings ----
# Signed test only (the undirected results stay in the csv/rds from 08-).
# One PDF per family of our rankings. Rows are their pooled gene sets
# (cell-type grouping x perturbation), columns are our rankings. Asterisks
# mark padj_family < 0.05. Crossed-out cells were not tested (fewer than
# minSize genes). Row bars: pooled DEG set sizes; column bars: our ranked
# genes; both split by one-to-one ortholog status.

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

nes_col <- colorRamp2(c(-3, 0, 3), c("#2a78d6", "#f0efec", "#e34948"))


# Load data ----
res_df <- readRDS(file.path(fld_in, "fgsea_lee_pooled_sets_vs_our_rankings.rds")) |>
  filter(version == "signed") |>
  mutate(
    sig = if_else(!is.na(padj_family) & padj_family < fdr_cut, "*", ""),
    row_id = paste(celltype, pert, sep = "|")
  )
ranking_meta <- readRDS(file.path(fld_in, "our_rankings.rds"))$meta
bd <- load_bar_data(fld_in)

pert_order <- sort(unique(res_df$pert))
row_df <- expand_grid(celltype = celltype_order, pert = pert_order) |>
  mutate(row_id = paste(celltype, pert, sep = "|"))


# Helper ----
plot_family <- function(rankings, title, file, width) {
  to_mat <- function(value_col, fill) {
    res_df |>
      filter(ranking %in% rankings) |>
      select(row_id, ranking, value = all_of(value_col)) |>
      complete(row_id = row_df$row_id, ranking = rankings, fill = list(value = fill)) |>
      pivot_wider(names_from = ranking, values_from = value) |>
      column_to_rownames("row_id") |>
      as.matrix() |>
      (\(m) m[row_df$row_id, rankings, drop = FALSE])()
  }
  nes <- to_mat("NES", NA_real_)
  tested <- to_mat("tested", FALSE)
  sig <- to_mat("sig", "")

  ht <- Heatmap(
    nes,
    name = "NES",
    col = nes_col,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(row_df$celltype, levels = celltype_order),
    row_labels = row_df$pert,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8),
    row_gap = unit(1.5, "mm"),
    column_labels = display_label(rankings),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 90,
    rect_gp = gpar(col = "#fcfcfb", lwd = 1),
    left_annotation = bar_anno(
      lee_set_counts(bd, row_df$celltype, row_df$pert, "all"),
      "row", "DEGs (up + down)"
    ),
    top_annotation = bar_anno(
      our_ranking_counts(bd, rankings),
      "column", "Ranked genes", size = unit(1.2, "cm")
    ),
    heatmap_legend_param = list(title = "NES", at = c(-3, -1.5, 0, 1.5, 3)),
    cell_fun = function(j, i, x, y, w, h, fill) {
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
  )

  pdf(file.path(fld_plot, file), width = width, height = 17)
  draw(
    ht,
    column_title = title,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    annotation_legend_list = list(
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
  )
  dev.off()
}


# Family heatmaps ----
plot_family(
  rankings = c("Domain-adjusted", "PRS"),
  title = "Pooled (up + down) DEG sets:\nlayer-adjusted SCZ and PRS",
  file = "heatmap_NES_pooled_domain_adjusted_PRS.pdf",
  width = 5.5
)

plot_family(
  rankings = ranking_meta$ranking[ranking_meta$family == "per_domain"],
  title = "Pooled (up + down) DEG sets: layer-restricted SCZ",
  file = "heatmap_NES_pooled_per_domain.pdf",
  width = 7
)

plot_family(
  rankings = ranking_meta$ranking[ranking_meta$family == "spg"],
  title = "Pooled (up + down) DEG sets: SPG-specific SCZ",
  file = "heatmap_NES_pooled_spg.pdf",
  width = 6
)


# Session info ----
sessioninfo::session_info()
