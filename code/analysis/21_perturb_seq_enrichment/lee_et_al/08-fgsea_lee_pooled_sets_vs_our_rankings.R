# fgsea: Lee et al. pooled (up + down) DEG sets against our rankings ----
# Same grid as 02-, but the up and down DEGs of each cell-type grouping x
# perturbation are pooled into one set ("<celltype>|<pert>|all"). Two tests:
# - signed: ranked by our t (scoreType = "std"). A set whose genes sit at both
#   ends of our ranking is diluted under this test.
# - undirected: ranked by |t| (scoreType = "pos"). Asks whether the DEGs are
#   concentrated among our most changed genes, whatever the direction.
# Universe as in 02-: genes tested in that mouse model x our ranked genes.
#
# Multiple testing: padj_ranking (BH within one ranking x version) and
# padj_family (BH within one family x version).

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
versions <- c(signed = "std", undirected = "pos")

fld_out <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")


# Load data ----
lee <- readRDS(file.path(fld_out, "lee_gene_sets.rds"))
universes <- lee$universes

ours <- readRDS(file.path(fld_out, "our_rankings.rds"))
our_rankings <- ours$rankings
ranking_meta <- ours$meta

ortho_121 <- read_csv(
  file.path(fld_out, "lee_ortholog_1to1.csv"),
  show_col_types = FALSE
)
ens2sym <- set_names(ortho_121$human_symbol, ortho_121$human_ensembl)


# Pooled gene sets ----
# Up and down sets are disjoint, so the union is their concatenation
set_universe_dir <- sub("\\|(up|down)$", "", names(lee$gene_sets))
pooled_sets <- split(lee$gene_sets, set_universe_dir) |>
  map(\(x) unique(unlist(x, use.names = FALSE)))
names(pooled_sets) <- paste0(names(pooled_sets), "|all")
stopifnot(sum(lengths(lee$gene_sets)) == sum(lengths(pooled_sets)))

set_universe <- sub("\\|all$", "", names(pooled_sets))
message(
  "Pooled sets: ", length(pooled_sets),
  "; above maxSize before universe matching: ", sum(lengths(pooled_sets) > max_size)
)


# Run fgsea ----
run_one <- function(ranking_name, universe_id, version) {
  score_type <- versions[[version]]
  stats <- our_rankings[[ranking_name]]
  stats <- stats[names(stats) %in% universes[[universe_id]]]
  if (version == "undirected") stats <- abs(stats)
  pathways <- pooled_sets[set_universe == universe_id]

  size_df <- tibble(
    set_id = names(pathways),
    n_set = lengths(pathways),
    n_in_ranking = map_int(pathways, ~ sum(.x %in% names(stats))),
    n_universe = length(stats)
  )

  res <- fgsea(
    pathways = pathways, stats = stats,
    minSize = min_size, maxSize = max_size, eps = 0,
    scoreType = score_type, BPPARAM = SerialParam()
  )
  res <- as_tibble(res) |> mutate(retried = FALSE)

  # Retry NA p-values (lopsided rankings) with more permutations
  failed <- res$pathway[is.na(res$pval)]
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

  size_df |>
    left_join(res |> rename(set_id = pathway), by = "set_id") |>
    mutate(version = version, ranking = ranking_name, .before = 1)
}

grid <- expand_grid(
  version = names(versions),
  ranking = names(our_rankings),
  universe_id = unique(set_universe)
)

set.seed(seed)
t0 <- Sys.time()
res_raw <- pmap(
  list(grid$ranking, grid$universe_id, grid$version),
  run_one,
  .progress = TRUE
) |>
  list_rbind()
message("fgsea elapsed: ", format(Sys.time() - t0))


# Tidy and adjust ----
bh_tested <- function(p, tested) {
  if_else(tested, p.adjust(if_else(tested, p, NA_real_), method = "BH"), NA_real_)
}

res_df <- res_raw |>
  separate_wider_delim(
    set_id,
    delim = "|",
    names = c("celltype", "pert", "dir"),
    cols_remove = FALSE
  ) |>
  left_join(ranking_meta |> select(ranking, family), by = "ranking") |>
  mutate(tested = !is.na(size)) |>
  group_by(version, ranking) |>
  mutate(padj_ranking = bh_tested(pval, tested)) |>
  group_by(version, family) |>
  mutate(padj_family = bh_tested(pval, tested)) |>
  ungroup() |>
  mutate(
    leading_edge_symbol = map_chr(
      leadingEdge,
      ~ if (is.null(.x)) NA_character_ else paste(ens2sym[.x], collapse = ";")
    )
  ) |>
  select(
    version, family, ranking, set_id, celltype, pert, dir,
    n_set, n_universe, n_in_ranking, size,
    ES, NES, pval, log2err, padj_ranking, padj_family,
    tested, retried, leadingEdge, leading_edge_symbol
  )


# Save ----
saveRDS(res_df, file.path(fld_out, "fgsea_lee_pooled_sets_vs_our_rankings.rds"))
write_csv(
  res_df |> select(-leadingEdge),
  file.path(fld_out, "fgsea_lee_pooled_sets_vs_our_rankings.csv")
)


# Summary ----
res_df |>
  group_by(version, family) |>
  summarise(
    n_tests = sum(tested),
    n_skipped_small = sum(!tested & n_in_ranking < min_size),
    n_skipped_large = sum(!tested & n_in_ranking > max_size),
    n_pval_na = sum(tested & is.na(pval)),
    n_fdr05 = sum(padj_family < 0.05, na.rm = TRUE),
    n_fdr05_neg = sum(padj_family < 0.05 & NES < 0, na.rm = TRUE),
    n_fdr05_pos = sum(padj_family < 0.05 & NES > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print()

res_df |>
  filter(tested) |>
  count(version, celltype, name = "n_tests") |>
  pivot_wider(names_from = version, values_from = n_tests) |>
  print()


# Session info ----
sessioninfo::session_info()
