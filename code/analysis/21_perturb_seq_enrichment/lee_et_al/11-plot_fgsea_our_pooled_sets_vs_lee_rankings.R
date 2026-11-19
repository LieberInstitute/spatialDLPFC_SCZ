# Plot fgsea results (secondary): our pooled DEG sets vs Lee et al. rankings ----
# Signed test only (the undirected results stay in the csv/rds from 10-).
# Rows are our pooled (up + down) DEG sets that passed the size filter in at
# least one test; columns are cell-type groupings x perturbations. Asterisks
# mark padj_family < 0.05. Crossed-out cells were not tested.
# Row bars: our pooled DEG set sizes; column bars: genes tested in their model
# (the ranking used); both split by one-to-one ortholog status.
#
# Usage: Rscript 11-plot_fgsea_our_pooled_sets_vs_lee_rankings.R [fdr10|p05]
# (same cutoff as the 10- run to plot; p05 reads the "_nomp05" outputs)

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

cutoff <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(cutoff)) cutoff <- "fdr10"
cutoff_cfg <- list(
  fdr10 = list(label = "FDR < 0.10", suffix = ""),
  p05 = list(label = "nominal p < 0.05", suffix = "_nomp05")
)[[cutoff]]
stopifnot(!is.null(cutoff_cfg))
sfx <- cutoff_cfg$suffix

fld_in <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")
fld_plot <- here("plots/21_perturb_seq_enrichment/lee_et_al")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)

celltype_order <- c(
  "Allcells", "Neuron",
  "Glutamatergic", "GABAergic", "Astro", "Oligo", "Microglia", "Vascular"
)
family_label <- c(
  domain_adjusted = "Layer-adjusted",
  per_domain = "Layer-restricted",
  spg = "SPG",
  prs = "PRS"
)

nes_col <- colorRamp2(c(-3, 0, 3), c("#2a78d6", "#f0efec", "#e34948"))


# Load data ----
res_df <- readRDS(file.path(fld_in, paste0("fgsea_our_pooled_sets_vs_lee_rankings", sfx, ".rds"))) |>
  filter(version == "signed") |>
  mutate(sig = !is.na(padj_family) & padj_family < fdr_cut)
ranking_meta <- readRDS(file.path(fld_in, "our_rankings.rds"))$meta
bd <- load_bar_data(fld_in)
# Our DEG gene lists at the chosen cutoff (written by 10-)
bd$our_degs <- read_csv(
  file.path(fld_in, paste0("our_pooled_deg_sets", sfx, ".csv")),
  show_col_types = FALSE
) |>
  mutate(matched = ensembl %in% bd$ortho_ens)

# Keep our sets tested at least once, in ranking order
row_df <- res_df |>
  filter(tested) |>
  distinct(family, ranking, set_id, n_set) |>
  mutate(ranking = factor(ranking, levels = ranking_meta$ranking)) |>
  arrange(ranking) |>
  mutate(row_label = sprintf("%s (n = %d)", display_label(as.character(ranking)), n_set))

col_df <- expand_grid(
  celltype = celltype_order,
  pert = sort(unique(res_df$pert))
) |>
  mutate(lee_id = paste(celltype, pert, sep = "|"))


# Matrices ----
to_mat <- function(value_col, fill) {
  res_df |>
    filter(set_id %in% row_df$set_id) |>
    select(set_id, lee_id, value = all_of(value_col)) |>
    complete(
      set_id = row_df$set_id,
      lee_id = col_df$lee_id,
      fill = list(value = fill)
    ) |>
    pivot_wider(names_from = lee_id, values_from = value) |>
    column_to_rownames("set_id") |>
    as.matrix() |>
    (\(m) m[row_df$set_id, col_df$lee_id, drop = FALSE])()
}

nes <- to_mat("NES", NA_real_)
tested <- to_mat("tested", FALSE)
sig <- to_mat("sig", FALSE)


# Heatmap ----
ht <- Heatmap(
  nes,
  name = "NES",
  col = nes_col,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(family_label[row_df$family], levels = family_label),
  row_labels = row_df$row_label,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 9),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),
  column_split = factor(col_df$celltype, levels = celltype_order),
  column_labels = col_df$pert,
  column_names_gp = gpar(fontsize = 7, fontface = "italic"),
  column_title_gp = gpar(fontsize = 9),
  column_gap = unit(2, "mm"),
  rect_gp = gpar(col = "#fcfcfb", lwd = 1),
  left_annotation = bar_anno(
    our_set_counts(bd, as.character(row_df$ranking), "all"),
    "row", "Our DEGs\n(up + down)", size = unit(2.6, "cm")
  ),
  top_annotation = bar_anno(
    lee_tested_counts(bd, col_df$celltype, col_df$pert),
    "column", "Tested genes", size = unit(1.2, "cm")
  ),
  heatmap_legend_param = list(title = "NES", at = c(-3, -1.5, 0, 1.5, 3)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!tested[i, j]) {
      grid.lines(
        x = unit.c(x - w / 2, x + w / 2),
        y = unit.c(y - h / 2, y + h / 2),
        gp = gpar(col = "#c3c2b7", lwd = 0.5)
      )
    } else if (sig[i, j]) {
      grid.text("*", x, y, vjust = 0.75, gp = gpar(fontsize = 9))
    }
  }
)

pdf(
  file.path(fld_plot, paste0("heatmap_NES_our_pooled_sets_vs_lee_rankings", sfx, ".pdf")),
  width = 17, height = 4 + 0.25 * nrow(row_df)
)
draw(
  ht,
  column_title = paste(
    paste0("Our pooled (up + down) SCZ-DEG sets (", cutoff_cfg$label, ") against Lee et al."),
    "perturbation rankings [sign(coef) x -log10(p)]"
  ),
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


# Session info ----
sessioninfo::session_info()
