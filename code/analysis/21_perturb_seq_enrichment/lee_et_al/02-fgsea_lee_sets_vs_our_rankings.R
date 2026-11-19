# fgsea: Lee et al. perturbation gene sets against our rankings ----
# For each of our rankings and each cell-type grouping x perturbation, the
# ranking is restricted to genes also tested in that mouse model (the
# universe), then the up and down DEG sets are tested with fgsea.
#
# Multiple testing:
# - padj_ranking: BH across all gene sets within one ranking
# - padj_family: BH across all rankings x gene sets within one family

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(fgsea)
  library(BiocParallel)
  library(sessioninfo)
})


# Config ----
min_size <- 15
max_size <- 1000
n_perm_retry <- 10000
seed <- 20261112

fld_out <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")


# Load data ----
lee <- readRDS(file.path(fld_out, "lee_gene_sets.rds"))
gene_sets <- lee$gene_sets
universes <- lee$universes

ours <- readRDS(file.path(fld_out, "our_rankings.rds"))
our_rankings <- ours$rankings
ranking_meta <- ours$meta

ortho_121 <- read_csv(
  file.path(fld_out, "lee_ortholog_1to1.csv"),
  show_col_types = FALSE
)
ens2sym <- set_names(ortho_121$human_symbol, ortho_121$human_ensembl)


# Run fgsea ----
set_universe <- sub("\\|(up|down)$", "", names(gene_sets))

run_one <- function(ranking_name, universe_id) {
  stats <- our_rankings[[ranking_name]]
  stats <- stats[names(stats) %in% universes[[universe_id]]]
  pathways <- gene_sets[set_universe == universe_id]

  # Record the size of every set after restricting to the ranked genes,
  # including those that fgsea will skip
  size_df <- tibble(
    set_id = names(pathways),
    n_set = lengths(pathways),
    n_in_ranking = map_int(pathways, ~ sum(.x %in% names(stats))),
    n_universe = length(stats)
  )

  res <- fgsea(
    pathways = pathways,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0,
    BPPARAM = SerialParam()
  )

  # When the ranking is lopsided in sign (e.g., SpD01-WMtz), fgsea can fail
  # to compute a p-value for large sets. Retry those with more permutations.
  failed <- res$pathway[is.na(res$pval)]
  if (length(failed) > 0) {
    res_retry <- fgsea(
      pathways = pathways[failed],
      stats = stats,
      minSize = min_size,
      maxSize = max_size,
      eps = 0,
      nPermSimple = n_perm_retry,
      BPPARAM = SerialParam()
    )
    res <- bind_rows(
      as_tibble(res) |> filter(!pathway %in% failed),
      as_tibble(res_retry)
    ) |>
      mutate(retried = pathway %in% failed)
  } else {
    res <- as_tibble(res) |> mutate(retried = FALSE)
  }

  size_df |>
    left_join(
      as_tibble(res) |> rename(set_id = pathway),
      by = "set_id"
    ) |>
    mutate(ranking = ranking_name, .before = 1)
}

grid <- expand_grid(
  ranking = names(our_rankings),
  universe_id = unique(set_universe)
)

set.seed(seed)
t0 <- Sys.time()
res_raw <- map2(
  grid$ranking, grid$universe_id,
  run_one,
  .progress = TRUE
) |>
  list_rbind()
message("fgsea elapsed: ", format(Sys.time() - t0))


# Tidy and adjust ----
res_df <- res_raw |>
  separate_wider_delim(
    set_id,
    delim = "|",
    names = c("celltype", "pert", "dir"),
    cols_remove = FALSE
  ) |>
  left_join(ranking_meta |> select(ranking, family), by = "ranking") |>
  # tested: the set passed the size filter (its p-value may still be NA)
  mutate(tested = !is.na(size)) |>
  group_by(ranking) |>
  mutate(padj_ranking = if_else(
    tested,
    p.adjust(if_else(tested, pval, NA_real_), method = "BH"),
    NA_real_
  )) |>
  group_by(family) |>
  mutate(padj_family = if_else(
    tested,
    p.adjust(if_else(tested, pval, NA_real_), method = "BH"),
    NA_real_
  )) |>
  ungroup() |>
  mutate(
    leading_edge_symbol = map_chr(
      leadingEdge,
      ~ if (is.null(.x)) NA_character_ else paste(ens2sym[.x], collapse = ";")
    ),
    leading_edge_ensembl = map_chr(
      leadingEdge,
      ~ if (is.null(.x)) NA_character_ else paste(.x, collapse = ";")
    )
  ) |>
  select(
    family, ranking, set_id, celltype, pert, dir,
    n_set, n_universe, n_in_ranking, size,
    ES, NES, pval, log2err, padj_fgsea = padj, padj_ranking, padj_family,
    tested, retried, leadingEdge, leading_edge_symbol, leading_edge_ensembl
  )

# fgsea's own padj only covers the two sets in one call, so drop it from
# the csv to avoid confusion
res_csv <- res_df |>
  select(-leadingEdge, -padj_fgsea)


# Save ----
saveRDS(res_df, file.path(fld_out, "fgsea_lee_sets_vs_our_rankings.rds"))
write_csv(res_csv, file.path(fld_out, "fgsea_lee_sets_vs_our_rankings.csv"))


# Summary ----
res_df |>
  group_by(family) |>
  summarise(
    n_tests = sum(tested),
    n_skipped_size = sum(!tested),
    n_pval_na = sum(tested & is.na(pval)),
    n_fdr05 = sum(padj_family < 0.05, na.rm = TRUE),
    n_fdr05_neg = sum(padj_family < 0.05 & NES < 0, na.rm = TRUE),
    n_fdr05_pos = sum(padj_family < 0.05 & NES > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print()

res_df |>
  filter(tested) |>
  count(celltype, dir, name = "n_tests") |>
  pivot_wider(names_from = dir, values_from = n_tests) |>
  print()


# Session info ----
sessioninfo::session_info()
