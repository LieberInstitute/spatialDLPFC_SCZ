# fgsea (secondary): our DEG sets against Lee et al. perturbation rankings ----
# Reverse direction of 02-. For each cell-type grouping x perturbation, their
# genes are ranked by sign(glm_coef) * -log10(glm_p). For each of our
# rankings, our up and down DEG sets (FDR < 0.10) are tested against their
# ranking, restricted to genes tested in both studies.
#
# Many of our sets are small and fall under minSize; treat this as a check.

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(fgsea)
  library(BiocParallel)
  library(sessioninfo)
})


# Config ----
our_fdr_cut <- 0.10
min_size <- 15
max_size <- 1000
n_perm_retry <- 10000
p_floor <- 1e-300 # glm_p can be exactly 0
seed <- 20261112

fld_out <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")


# Load data ----
glm_human <- readRDS(file.path(fld_out, "lee_glm_human.rds"))
ranking_meta <- readRDS(file.path(fld_out, "our_rankings.rds"))$meta
our_rankings <- readRDS(file.path(fld_out, "our_rankings.rds"))$rankings

spd_order <- ranking_meta$ranking[ranking_meta$family == "per_domain"]


# Our DEG sets ----
# Need the FDR column for each ranking, so reload the sources here.
## Domain-adjusted ----
domain_adjusted_df <- read_csv(
  here("processed-data/rds/10_dx_deg_adjust_spd", "dx-deg_PRECAST07.csv"),
  show_col_types = FALSE
) |>
  transmute(ranking = "Domain-adjusted", ensembl, t = t_stat_scz, fdr = fdr_scz)

## Per-domain ----
per_domain_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  ),
  show_col_types = FALSE
) |>
  transmute(ranking = PRECAST_spd, ensembl = gene_id, t, fdr = adj.P.Val)

## SPG ----
spg_files <- c(
  NeuN = "dx_deg_spg_neun/neun-dx_DEG-GM.csv",
  Neuropil = "dx_deg_spg_neuropil/neuropil-dx_DEG-GM.csv",
  PNN = "dx_deg_spg_pnn/pnn-dx_DEG-GM.csv",
  Vasc = "dx_deg_spg_vasc/vasc-dx_DEG-allSpD.csv"
)
spg_df <- imap(spg_files, ~ read_csv(
  here("processed-data/20_dx_deg_spg", .x),
  show_col_types = FALSE
) |>
  transmute(ranking = .y, ensembl, t = t_stat_scz, fdr = fdr_scz)) |>
  list_rbind()

## PRS ----
prs_df <- read_csv(
  here(
    "processed-data/rds/14_prs_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  ),
  show_col_types = FALSE
) |>
  transmute(ranking = "PRS", ensembl = gene_id, t, fdr = adj.P.Val)

our_df <- bind_rows(domain_adjusted_df, per_domain_df, spg_df, prs_df) |>
  filter(!is.na(fdr), !is.na(t))
stopifnot(setequal(unique(our_df$ranking), ranking_meta$ranking))

our_deg_df <- our_df |>
  filter(fdr < our_fdr_cut) |>
  mutate(
    dir = if_else(t > 0, "up", "down"),
    set_id = paste(ranking, dir, sep = "|")
  )

our_sets <- split(our_deg_df$ensembl, our_deg_df$set_id)

our_set_size <- our_deg_df |>
  count(ranking, dir, name = "n_deg") |>
  complete(ranking = ranking_meta$ranking, dir = c("down", "up"), fill = list(n_deg = 0))
print(our_set_size, n = Inf)


# Their rankings ----
lee_rank_df <- glm_human |>
  mutate(
    score = sign(glm_coef) * -log10(pmax(glm_p, p_floor)),
    lee_id = paste(celltype, pert, sep = "|")
  )

lee_rankings <- split(
  set_names(lee_rank_df$score, lee_rank_df$human_ensembl),
  lee_rank_df$lee_id
)


# Run fgsea ----
run_one <- function(our_ranking, lee_id) {
  # Universe: genes tested in both our ranking and their model
  stats <- lee_rankings[[lee_id]]
  stats <- stats[names(stats) %in% names(our_rankings[[our_ranking]])]

  pathways <- our_sets[startsWith(names(our_sets), paste0(our_ranking, "|"))]
  if (length(pathways) == 0) {
    return(NULL)
  }

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
    left_join(res |> rename(set_id = pathway), by = "set_id") |>
    mutate(lee_id = lee_id, .before = 1)
}

grid <- expand_grid(
  our_ranking = ranking_meta$ranking,
  lee_id = names(lee_rankings)
)

set.seed(seed)
t0 <- Sys.time()
res_raw <- map2(
  grid$our_ranking, grid$lee_id,
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
    names = c("ranking", "dir"),
    cols_remove = FALSE
  ) |>
  separate_wider_delim(
    lee_id,
    delim = "|",
    names = c("celltype", "pert"),
    cols_remove = FALSE
  ) |>
  left_join(ranking_meta |> select(ranking, family), by = "ranking") |>
  mutate(tested = !is.na(size)) |>
  group_by(family) |>
  mutate(padj_family = if_else(
    tested,
    p.adjust(if_else(tested, pval, NA_real_), method = "BH"),
    NA_real_
  )) |>
  ungroup() |>
  mutate(
    leading_edge_ensembl = map_chr(
      leadingEdge,
      ~ if (is.null(.x)) NA_character_ else paste(.x, collapse = ";")
    )
  ) |>
  select(
    family, ranking, dir, set_id, celltype, pert, lee_id,
    n_set, n_universe, n_in_ranking, size,
    ES, NES, pval, log2err, padj_family,
    tested, retried, leadingEdge, leading_edge_ensembl
  )


# Save ----
write_csv(our_set_size, file.path(fld_out, "our_deg_set_size.csv"))
# Gene lists of our DEG sets (for the heatmap bar annotations)
write_csv(
  our_deg_df |> select(ranking, dir, ensembl, t, fdr),
  file.path(fld_out, "our_deg_sets.csv")
)
saveRDS(res_df, file.path(fld_out, "fgsea_our_sets_vs_lee_rankings.rds"))
write_csv(
  res_df |> select(-leadingEdge),
  file.path(fld_out, "fgsea_our_sets_vs_lee_rankings.csv")
)


# Summary ----
res_df |>
  group_by(family, ranking, dir) |>
  summarise(
    n_set = first(n_set),
    n_tests = sum(tested),
    n_fdr05 = sum(padj_family < 0.05, na.rm = TRUE),
    n_fdr05_neg = sum(padj_family < 0.05 & NES < 0, na.rm = TRUE),
    n_fdr05_pos = sum(padj_family < 0.05 & NES > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print(n = Inf)


# Session info ----
sessioninfo::session_info()
