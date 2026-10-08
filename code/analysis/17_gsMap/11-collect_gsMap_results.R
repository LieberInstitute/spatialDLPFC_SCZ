# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
fld_gsmap <- here("processed-data/17_gsMap")
fld_workdir <- file.path(fld_gsmap, "workdir")
fld_out <- here("processed-data/rds/17_gsMap")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)

sample_ids <- readLines(file.path(fld_gsmap, "h5ad", "sample_list.txt"))

trait_names <- read_tsv(
  here("code/analysis/17_gsMap/traits.tsv"),
  comment = "#",
  show_col_types = FALSE
) |>
  pull(trait_name)

cat("Traits:", trait_names, "\n")


# Load data ----
## Spot-level metadata from the finalized spe ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

spot_meta <- colData(spe) |>
  as.data.frame() |>
  select(key, sample_id, brnum, dx, age, sex, spd_label) |>
  as_tibble()

rm(spe)


# Spot-level spatial LDSC results ----
## {workdir}/{sample}/spatial_ldsc/{sample}_{trait}.csv.gz; `spot` = spe$key
read_spot_ldsc <- function(.smp, .trait) {
  file.path(
    fld_workdir, .smp, "spatial_ldsc",
    paste0(.smp, "_", .trait, ".csv.gz")
  ) |>
    read_csv(show_col_types = FALSE) |>
    select(key = spot, any_of("beta"), p) |>
    mutate(key = as.character(key), trait = .trait, .after = key)
}

spot_res_df <- expand_grid(sample_id = sample_ids, trait = trait_names) |>
  pmap(.f = function(sample_id, trait) read_spot_ldsc(sample_id, trait)) |>
  list_rbind() |>
  mutate(neg_log10_p = -log10(p)) |>
  inner_join(spot_meta, by = "key", relationship = "many-to-one")

# error prevention: every spot has a result for every trait
stopifnot(nrow(spot_res_df) == nrow(spot_meta) * length(trait_names))
stopifnot(!anyNA(spot_res_df$p))


# Sample x domain summaries ----
## gsMap per-sample Cauchy combination (08) ----
read_sample_cauchy <- function(.smp, .trait) {
  file.path(
    fld_workdir, .smp, "cauchy_combination",
    paste0(.smp, "_", .trait, ".Cauchy.csv.gz")
  ) |>
    read_csv(show_col_types = FALSE) |>
    transmute(
      sample_id = .smp,
      trait = .trait,
      spd_label = annotation,
      p_cauchy,
      p_median
    )
}

sample_cauchy_df <- expand_grid(sample_id = sample_ids, trait = trait_names) |>
  pmap(.f = function(sample_id, trait) read_sample_cauchy(sample_id, trait)) |>
  list_rbind()

## Spot-level summaries per sample x domain ----
## Mean -log10(p) is less sensitive to the number of spots in a domain than
## the Cauchy p-value, so both are kept for donor-level comparisons.
sample_spd_df <- spot_res_df |>
  group_by(trait, sample_id, brnum, dx, age, sex, spd_label) |>
  summarize(
    n_spots = n(),
    mean_neg_log10_p = mean(neg_log10_p),
    prop_p_lt_0.05 = mean(p < 0.05),
    .groups = "drop"
  ) |>
  mutate(spd_label = as.character(spd_label)) |>
  left_join(
    sample_cauchy_df,
    by = c("trait", "sample_id", "spd_label"),
    relationship = "one-to-one"
  ) |>
  mutate(neg_log10_p_cauchy = -log10(p_cauchy))

# error prevention
stopifnot(!anyNA(sample_spd_df$p_cauchy))


# Across-sample Cauchy combination (09) ----
across_cauchy_df <- expand_grid(trait = trait_names, group = c("all", "ntc", "scz")) |>
  pmap(.f = function(trait, group) {
    file.path(
      fld_gsmap, "cauchy_across_samples",
      paste0(trait, "_", group, ".Cauchy.csv.gz")
    ) |>
      read_csv(show_col_types = FALSE) |>
      transmute(trait = trait, group = group, spd_label = annotation, p_cauchy, p_median)
  }) |>
  list_rbind() |>
  mutate(neg_log10_p_cauchy = -log10(p_cauchy))


# Save ----
# NOTE: spot-level results are kept slim and key-indexed so they can be
# joined back into colData(spe), following 16_integration_ling.
saveRDS(
  spot_res_df |> select(key, sample_id, trait, any_of("beta"), p, neg_log10_p),
  file.path(fld_out, "gsMap_spot_res_df.rds")
)
saveRDS(sample_spd_df, file.path(fld_out, "gsMap_sample_spd_df.rds"))
saveRDS(across_cauchy_df, file.path(fld_out, "gsMap_across_sample_cauchy_df.rds"))

write_csv(sample_spd_df, file.path(fld_out, "gsMap_sample_spd_summary.csv"))
write_csv(across_cauchy_df, file.path(fld_out, "gsMap_across_sample_cauchy.csv"))

print("Finished collecting gsMap results")


# Session info ----
sessioninfo::session_info()
