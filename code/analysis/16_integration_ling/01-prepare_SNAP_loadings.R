# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
path_snap_xlsx <- here(
  "code/analysis/16_integration_ling",
  "Ling2024_SNAP_ready_to_use.xlsx"
)

fld_out <- here("processed-data/rds/16_integration_ling")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)


# Load data ----
## Full, descending-ranked cNMF gene loadings ----
snapa_full <- read_excel(path_snap_xlsx, sheet = "SNAPa_loadings_ranked") |>
  rename(loading = `SNAP-a`) |>
  mutate(program = "SNAP-a", .after = gene)

snapn_full <- read_excel(path_snap_xlsx, sheet = "SNAPn_loadings_ranked") |>
  rename(loading = `SNAP-n`) |>
  mutate(program = "SNAP-n", .after = gene)

## Top 2000 genes by loading - Ling et al.'s scoring convention ----
snapa_top2000 <- read_excel(path_snap_xlsx, sheet = "SNAPa_top2000") |>
  rename(loading = `SNAP-a`) |>
  mutate(program = "SNAP-a", .after = gene)

snapn_top2000 <- read_excel(path_snap_xlsx, sheet = "SNAPn_top2000") |>
  rename(loading = `SNAP-n`) |>
  mutate(program = "SNAP-n", .after = gene)

# NOTE: the xlsx also contains `donor_SNAP_scores` and `donor_metadata`
# sheets, but these are Ling et al.'s own published per-donor scores/
# metadata from their cohort (not ours), and are not used anywhere in this
# projection-onto-Visium analysis, so they are intentionally not loaded here.


# Sanity checks ----
stopifnot(
  "cNMF loadings are expected to be non-negative" =
    all(snapa_top2000$loading >= 0) &&
      all(snapn_top2000$loading >= 0)
)

stopifnot(
  "Duplicated gene symbols found in top 2000 gene set" =
    !any(duplicated(snapa_top2000$gene)) &&
      !any(duplicated(snapn_top2000$gene))
)

stopifnot(nrow(snapa_top2000) == 2000)
stopifnot(nrow(snapn_top2000) == 2000)


# Combine and save ----
snap_loadings <- list(
  full = bind_rows(snapa_full, snapn_full),
  top2000 = bind_rows(snapa_top2000, snapn_top2000)
)

saveRDS(
  snap_loadings,
  file.path(fld_out, "SNAP_loadings.rds")
)

print("Finished preparing Ling et al. SNAP loadings")


# Session info ----
sessioninfo::session_info()
