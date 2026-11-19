# Fig 6C check with all DEGs (up and down pooled) ----
# Allcells Matr3 and Scyl1 DEGs (fdr < 0.05, one-to-one orthologs), with up
# and down pooled into one set, tested against PRS and per-domain rankings.
# Two tests per set:
# - Signed: ranked by our t, as in 02. A pooled set whose genes move in both
#   directions can partly cancel out under this test.
# - Undirected: ranked by |t| with scoreType = "pos". Asks whether the DEGs
#   are concentrated among our most changed genes, whatever the direction.
# Universe as in 02: genes tested in that mouse model x our ranked genes.
# BH is applied within each test type across the 2 sets x 8 rankings.

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(fgsea)
  library(BiocParallel)
  library(ComplexHeatmap)
  library(circlize)
  library(sessioninfo)
})


# Config ----
min_size <- 15
max_size <- 1000
n_perm_retry <- 10000
seed <- 20261112
fdr_cut <- 0.05
perts <- c("Matr3", "Scyl1")

fld_in <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")
fld_plot <- here("plots/21_perturb_seq_enrichment/lee_et_al")


# Load data ----
lee <- readRDS(file.path(fld_in, "lee_gene_sets.rds"))
ours <- readRDS(file.path(fld_in, "our_rankings.rds"))
ranking_meta <- ours$meta

rankings <- c("PRS", ranking_meta$ranking[ranking_meta$family == "per_domain"])

# Pooled sets: union of up and down (they are disjoint by construction)
pooled_sets <- map(perts, \(p) {
  union(
    lee$gene_sets[[paste("Allcells", p, "down", sep = "|")]],
    lee$gene_sets[[paste("Allcells", p, "up", sep = "|")]]
  )
}) |>
  set_names(perts)
print(lengths(pooled_sets))


# Run fgsea ----
run_fgsea <- function(pathways, stats, score_type) {
  res <- fgsea(
    pathways = pathways, stats = stats,
    minSize = min_size, maxSize = max_size, eps = 0,
    scoreType = score_type, BPPARAM = SerialParam()
  )
  failed <- res$pathway[is.na(res$pval)]
  res <- as_tibble(res) |> mutate(retried = FALSE)
  if (length(failed) > 0) {
    res_retry <- fgsea(
      pathways = pathways[failed], stats = stats,
      minSize = min_size, maxSize = max_size, eps = 0,
      scoreType = score_type, nPermSimple = n_perm_retry,
      BPPARAM = SerialParam()
    )
    res <- bind_rows(
      res |> filter(!pathway %in% failed),
      as_tibble(res_retry) |> mutate(retried = TRUE)
    )
  }
  res
}

grid <- expand_grid(
  test = c("Signed (t)", "Undirected (|t|)"),
  ranking = rankings,
  pert = perts
)

set.seed(seed)
res_df <- pmap(grid, \(test, ranking, pert) {
  stats <- ours$rankings[[ranking]]
  stats <- stats[names(stats) %in% lee$universes[[paste("Allcells", pert, sep = "|")]]]
  score_type <- "std"
  if (test == "Undirected (|t|)") {
    stats <- abs(stats)
    score_type <- "pos"
  }
  run_fgsea(pooled_sets[pert], stats, score_type) |>
    mutate(test = test, ranking = ranking, n_universe = length(stats))
}) |>
  list_rbind() |>
  rename(pert = pathway) |>
  group_by(test) |>
  mutate(padj = p.adjust(pval, method = "BH")) |>
  ungroup() |>
  mutate(sig = if_else(!is.na(padj) & padj < fdr_cut, "*", ""))

res_df |>
  select(test, pert, ranking, size, NES, pval, padj, retried) |>
  arrange(test, pert, match(ranking, rankings)) |>
  print(n = Inf)

write_csv(
  res_df |> select(test, pert, ranking, n_universe, size, ES, NES, pval, padj, retried),
  file.path(fld_in, "fig6C_check_all_DEGs.csv")
)


# Plot ----
# Signed test only (both tests are in the csv). Row bars: pooled DEG set
# sizes; column bars: our ranked genes; both split by one-to-one ortholog
# status.
source(here("code/analysis/21_perturb_seq_enrichment/lee_et_al/utils_heatmap_bars.R"))
bd <- load_bar_data(fld_in)

nes_col <- colorRamp2(c(-3, 0, 3), c("#2a78d6", "#f0efec", "#e34948"))

plot_df <- res_df |>
  filter(test == "Signed (t)") |>
  mutate(ranking = factor(ranking, levels = rankings))
to_mat <- function(v) {
  plot_df |>
    select(pert, ranking, value = all_of(v)) |>
    pivot_wider(names_from = ranking, values_from = value) |>
    column_to_rownames("pert") |>
    as.matrix() |>
    (\(m) m[perts, rankings, drop = FALSE])()
}
nes <- to_mat("NES")
sig <- to_mat("sig")

ht <- Heatmap(
  nes,
  name = "NES",
  col = nes_col,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(
    if_else(rankings == "PRS", "PRS", "Layer-restricted SCZ"),
    levels = c("PRS", "Layer-restricted SCZ")
  ),
  column_gap = unit(2, "mm"),
  row_labels = paste(perts, "up + down"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "#fcfcfb", lwd = 1),
  left_annotation = bar_anno(
    lee_set_counts(bd, "Allcells", perts, "all"),
    "row", "DEGs"
  ),
  top_annotation = bar_anno(
    our_ranking_counts(bd, rankings),
    "column", "Ranked genes", size = unit(1.2, "cm")
  ),
  heatmap_legend_param = list(title = "NES", at = c(-3, -1.5, 0, 1.5, 3)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    label <- if (is.na(nes[i, j])) "NA" else sprintf("%.1f%s", nes[i, j], sig[i, j])
    grid.text(label, x, y, gp = gpar(fontsize = 7))
  }
)

pdf(file.path(fld_plot, "fig6C_check_Allcells_Matr3_Scyl1_all_DEGs.pdf"), width = 7.5, height = 3)
draw(
  ht,
  column_title = "Allcells Matr3 / Scyl1 DEGs, up and down pooled (cf. Lee et al. Fig 6C)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_list = list(bar_legend())
)
dev.off()


# Session info ----
sessioninfo::session_info()
