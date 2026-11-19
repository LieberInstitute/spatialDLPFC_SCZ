# Helpers: stacked bar annotations of gene-set sizes for the NES heatmaps ----
# Each bar is split into genes with a one-to-one ortholog (the mapping used in
# the analysis, lee_ortholog_1to1.csv) and genes without one.
# - Their genes (mouse): one-to-one human ortholog or not.
# - Our genes (human): one-to-one match to a mouse gene tested in CROP-seq
#   or not.
# Sourced by the plotting scripts (03, 05, 07, 09, 11). Expects the packages
# loaded by those scripts (tidyverse, ComplexHeatmap).

# Display names for our rankings in the figures. Data files keep the
# original names ("Domain-adjusted", family ids domain_adjusted/per_domain).
display_label <- function(x) {
  str_replace_all(x, c(
    "Domain-adjusted" = "Layer-adjusted",
    "domain-adjusted" = "layer-adjusted",
    "Per-domain" = "Layer-restricted",
    "per-domain" = "layer-restricted"
  ))
}

bar_cols <- c(matched = "#5f5e57", unmatched = "#d6d5cc")
bar_labels <- c(matched = "1-1 ortholog", unmatched = "No 1-1 ortholog")


# Data ----
load_bar_data <- function(fld_in) {
  ortho_ens <- read_csv(
    file.path(fld_in, "lee_ortholog_1to1.csv"),
    show_col_types = FALSE
  )$human_ensembl

  list(
    # Their DEGs: one row per DEG, human_ensembl is NA without a 1-1 ortholog
    lee_deg = read_csv(
      file.path(fld_in, "lee_deg_mapping.csv"),
      show_col_types = FALSE
    ) |>
      mutate(matched = !is.na(human_ensembl)),
    # Their tested genes per model
    lee_tested = read_csv(
      file.path(fld_in, "lee_gene_set_size.csv"),
      show_col_types = FALSE
    ) |>
      distinct(celltype, pert, n_tested_mouse, n_tested_human_1to1),
    our_rankings = readRDS(file.path(fld_in, "our_rankings.rds"))$rankings,
    our_degs = read_csv(
      file.path(fld_in, "our_deg_sets.csv"),
      show_col_types = FALSE
    ) |>
      mutate(matched = ensembl %in% ortho_ens),
    ortho_ens = ortho_ens
  )
}


# Count matrices (rows = sets, columns = matched, unmatched) ----
to_count_mat <- function(matched, unmatched, ids) {
  matrix(
    c(matched, unmatched),
    ncol = 2,
    dimnames = list(ids, names(bar_cols))
  )
}

# Their DEG sets. dir = "all" pools up and down. Sets without DEGs get 0.
lee_set_counts <- function(bd, celltype, pert, dir) {
  ids <- paste(celltype, pert, dir, sep = "|")
  counts <- bind_rows(bd$lee_deg, bd$lee_deg |> mutate(dir = "all")) |>
    group_by(id = paste(celltype, pert, dir, sep = "|")) |>
    # unmatched first: summarise() reuses a column name once it is overwritten
    summarise(unmatched = sum(!matched), matched = sum(matched), .groups = "drop")
  counts <- counts[match(ids, counts$id), ]
  to_count_mat(
    replace_na(counts$matched, 0L),
    replace_na(counts$unmatched, 0L),
    ids
  )
}

# Their tested genes per model (the rankings used in the reverse analysis)
lee_tested_counts <- function(bd, celltype, pert) {
  counts <- tibble(celltype, pert) |>
    left_join(bd$lee_tested, by = c("celltype", "pert"))
  stopifnot(!anyNA(counts$n_tested_mouse))
  to_count_mat(
    counts$n_tested_human_1to1,
    counts$n_tested_mouse - counts$n_tested_human_1to1,
    paste(celltype, pert, sep = "|")
  )
}

# Our rankings: all ranked genes
our_ranking_counts <- function(bd, rankings) {
  matched <- map_int(rankings, ~ sum(names(bd$our_rankings[[.x]]) %in% bd$ortho_ens))
  total <- map_int(rankings, ~ length(bd$our_rankings[[.x]]))
  to_count_mat(matched, total - matched, rankings)
}

# Our DEG sets. dir = "all" pools up and down.
our_set_counts <- function(bd, ranking, dir) {
  counts <- map2(ranking, dir, \(r, d) {
    g <- bd$our_degs |> filter(ranking == r, d == "all" | dir == d)
    c(sum(g$matched), sum(!g$matched))
  })
  to_count_mat(
    map_dbl(counts, 1), map_dbl(counts, 2),
    paste(ranking, dir, sep = "|")
  )
}


# Annotations ----
# Row bars grow leftward from the heatmap (0 next to the cells), with the
# name below. Column bars have the name on the left and the axis on the
# right, so neither collides with the other annotation.
bar_anno <- function(mat, which = c("row", "column"), title, size = unit(1.8, "cm")) {
  which <- match.arg(which)
  if (which == "row") {
    rowAnnotation(
      n = anno_barplot(
        mat,
        which = "row",
        gp = gpar(fill = bar_cols, col = NA),
        bar_width = 0.8,
        border = FALSE,
        axis_param = list(gp = gpar(fontsize = 6), direction = "reverse"),
        width = size
      ),
      annotation_label = title,
      annotation_name_gp = gpar(fontsize = 7),
      annotation_name_rot = 0,
      annotation_name_side = "bottom",
      annotation_name_offset = unit(6, "mm")
    )
  } else {
    HeatmapAnnotation(
      n = anno_barplot(
        mat,
        which = "column",
        gp = gpar(fill = bar_cols, col = NA),
        bar_width = 0.8,
        border = FALSE,
        axis_param = list(gp = gpar(fontsize = 6), side = "right"),
        height = size
      ),
      annotation_label = title,
      annotation_name_gp = gpar(fontsize = 7),
      annotation_name_side = "left"
    )
  }
}

bar_legend <- function() {
  Legend(
    labels = bar_labels,
    title = "Genes",
    legend_gp = gpar(fill = bar_cols)
  )
}
