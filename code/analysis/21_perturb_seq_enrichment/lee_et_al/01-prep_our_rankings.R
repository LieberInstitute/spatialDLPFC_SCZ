# Collect our SCZ-DEG and PRS rankings ----
# Each ranking is a named numeric vector of moderated t-statistics, named by
# human Ensembl ID. Families:
# - domain_adjusted: SCZ vs NTC, adjusted for spatial domain (1 ranking)
# - per_domain: SCZ vs NTC within each spatial domain (7 rankings)
# - spg: SCZ vs NTC within SPG-positive spots (4 rankings)
# - prs: association with normalized PRS (1 ranking)

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
fld_out <- here("processed-data/21_perturb_seq_enrichment/lee_et_al")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)

spd_order <- c(
  "SpD07-L1/M",
  "SpD06-L2/3",
  "SpD02-L3/4",
  "SpD05-L5",
  "SpD03-L6",
  "SpD01-WMtz",
  "SpD04-WM"
)


# Helper ----
# Turn a data frame into a named t-statistic vector, with checks
to_ranking <- function(df, id_col, stat_col) {
  df <- df |>
    filter(!is.na(.data[[id_col]]), !is.na(.data[[stat_col]]))
  stopifnot(!anyDuplicated(df[[id_col]]))
  set_names(df[[stat_col]], df[[id_col]])
}


# Load rankings ----
## Domain-adjusted SCZ-DEGs ----
domain_adjusted <- read_csv(
  here("processed-data/rds/10_dx_deg_adjust_spd", "dx-deg_PRECAST07.csv"),
  show_col_types = FALSE
) |>
  to_ranking("ensembl", "t_stat_scz")

## Per-domain SCZ-DEGs ----
per_domain_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  ),
  show_col_types = FALSE
)
stopifnot(setequal(unique(per_domain_df$PRECAST_spd), spd_order))

per_domain <- spd_order |>
  set_names() |>
  map(~ per_domain_df |>
    filter(PRECAST_spd == .x) |>
    to_ranking("gene_id", "t"))

## SPG-specific SCZ-DEGs ----
spg_files <- c(
  NeuN = "dx_deg_spg_neun/neun-dx_DEG-GM.csv",
  Neuropil = "dx_deg_spg_neuropil/neuropil-dx_DEG-GM.csv",
  PNN = "dx_deg_spg_pnn/pnn-dx_DEG-GM.csv",
  Vasc = "dx_deg_spg_vasc/vasc-dx_DEG-allSpD.csv"
)

spg <- spg_files |>
  map(~ read_csv(
    here("processed-data/20_dx_deg_spg", .x),
    show_col_types = FALSE
  ) |>
    to_ranking("ensembl", "t_stat_scz"))

## PRS-associated genes ----
prs <- read_csv(
  here(
    "processed-data/rds/14_prs_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  ),
  show_col_types = FALSE
) |>
  to_ranking("gene_id", "t")


# Combine ----
our_rankings <- c(
  list(`Domain-adjusted` = domain_adjusted),
  per_domain,
  spg,
  list(PRS = prs)
)

ranking_meta <- tibble(
  ranking = names(our_rankings),
  family = c(
    "domain_adjusted",
    rep("per_domain", length(per_domain)),
    rep("spg", length(spg)),
    "prs"
  ),
  n_genes = lengths(our_rankings),
  n_pos = map_int(our_rankings, ~ sum(.x > 0)),
  n_neg = map_int(our_rankings, ~ sum(.x < 0))
)
print(ranking_meta, n = Inf)

saveRDS(
  list(rankings = our_rankings, meta = ranking_meta),
  file.path(fld_out, "our_rankings.rds")
)


# Session info ----
sessioninfo::session_info()
