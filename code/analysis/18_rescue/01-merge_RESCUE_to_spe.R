# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
fld_out <- here("processed-data/rds/18_rescue")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)

# Non-proportion (metadata) columns of the RESCUE spot-level output; every
# other column is a cell type proportion.
rescue_meta_cols <- c(
  "spot_id", "sample_id", "barcode", "array_row", "array_col", "domain",
  "pxl_col_in_fullres", "pxl_row_in_fullres", "nUMI", "unattributed_ratio"
)

# Tolerance for proportions summing to 1 within a spot
prop_sum_tol <- 1e-6


# Load data ----
## Finalized spe (63 donors, kept/QC-passed spots only) ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

## RESCUE spot-level deconvolution results ----
rescue_df <- read_csv(
  here("processed-data/rds/RESCUE_res", "spot_level_results.csv"),
  show_col_types = FALSE
)

ct_cols <- setdiff(colnames(rescue_df), rescue_meta_cols)

cat(sprintf(
  "RESCUE: %d rows, %d cell type columns (incl. 'Unattributed Expression')\n",
  nrow(rescue_df), length(ct_cols)
))
print(ct_cols)


# Merge QC -------------------------------------------------------------------
## Each check is recorded as a row of `qc_df`; `hard = TRUE` checks indicate
## misalignment and stop the script if they fail. `hard = FALSE` checks are
## documented diagnostics only.
qc_list <- list()
add_qc <- function(check, value, pass, hard, note = "") {
  qc_list[[length(qc_list) + 1]] <<- tibble(
    check = check, value = as.character(value),
    pass = pass, hard = hard, note = note
  )
}

## 1) Key uniqueness & coverage ----
# `spot_id` in RESCUE is "<barcode>_<sample_id>", identical to `spe$key`
stopifnot(identical(colnames(spe), spe$key))

n_spe <- ncol(spe)
n_rescue <- nrow(rescue_df)
n_dup_rescue <- sum(duplicated(rescue_df$spot_id))
n_rescue_in_spe <- sum(rescue_df$spot_id %in% spe$key)
n_spe_in_rescue <- sum(spe$key %in% rescue_df$spot_id)

add_qc("n_spots_spe", n_spe, TRUE, FALSE)
add_qc("n_spots_rescue", n_rescue, TRUE, FALSE)
add_qc("n_duplicated_spot_id_rescue", n_dup_rescue, n_dup_rescue == 0, TRUE)
add_qc(
  "n_rescue_spots_matched_to_spe", n_rescue_in_spe,
  n_rescue_in_spe == n_rescue, FALSE,
  sprintf("%.4f%% of RESCUE spots", 100 * n_rescue_in_spe / n_rescue)
)
add_qc(
  "n_spe_spots_with_rescue", n_spe_in_rescue,
  n_spe_in_rescue == n_spe, FALSE,
  sprintf("%.4f%% of spe spots (unmatched spots get NA)", 100 * n_spe_in_rescue / n_spe)
)
add_qc(
  "identical_spot_id_and_barcode_cols",
  identical(rescue_df$spot_id, rescue_df$barcode),
  identical(rescue_df$spot_id, rescue_df$barcode), FALSE
)

## 2) Sample coverage ----
smp_spe <- unique(spe$sample_id)
smp_rescue <- unique(rescue_df$sample_id)
add_qc("n_samples_spe", length(smp_spe), TRUE, FALSE)
add_qc("n_samples_rescue", length(smp_rescue), TRUE, FALSE)
add_qc(
  "samples_only_in_spe", paste(setdiff(smp_spe, smp_rescue), collapse = ";"),
  length(setdiff(smp_spe, smp_rescue)) == 0, FALSE
)
add_qc(
  "samples_only_in_rescue", paste(setdiff(smp_rescue, smp_spe), collapse = ";"),
  length(setdiff(smp_rescue, smp_spe)) == 0, FALSE
)

## 3) Concordance of shared metadata among matched spots ----
idx <- match(rescue_df$spot_id, spe$key)
ok <- !is.na(idx)
rescue_ok <- rescue_df[ok, ]
spe_cd <- colData(spe)[idx[ok], ]
spe_xy <- spatialCoords(spe)[idx[ok], , drop = FALSE]

concord_frac <- function(x, y) mean(as.character(x) == as.character(y))

concord <- c(
  sample_id = concord_frac(rescue_ok$sample_id, spe_cd$sample_id),
  array_row = concord_frac(rescue_ok$array_row, spe_cd$array_row),
  array_col = concord_frac(rescue_ok$array_col, spe_cd$array_col),
  pxl_row_in_fullres = concord_frac(
    rescue_ok$pxl_row_in_fullres, spe_xy[, "pxl_row_in_fullres"]
  ),
  pxl_col_in_fullres = concord_frac(
    rescue_ok$pxl_col_in_fullres, spe_xy[, "pxl_col_in_fullres"]
  ),
  # RESCUE `domain` (spd01-spd07) is the PRECAST k = 7 label
  domain_vs_PRECAST_07 = concord_frac(rescue_ok$domain, spe_cd$PRECAST_07)
)

iwalk(concord, \(frac, var) {
  add_qc(
    sprintf("concordance_%s", var), sprintf("%.6f", frac),
    frac == 1, TRUE, "fraction of matched spots with identical value"
  )
})

## 4) nUMI vs spe$sum_umi (diagnostic only) ----
# RESCUE's nUMI is computed on the genes it used, so it is expected to be
# slightly smaller than, but highly correlated with, total UMI in `spe`.
umi_ratio <- rescue_ok$nUMI / spe_cd$sum_umi
add_qc("nUMI_exact_match_frac", sprintf("%.4f", mean(rescue_ok$nUMI == spe_cd$sum_umi)), NA, FALSE)
add_qc("nUMI_pearson_cor_with_sum_umi", sprintf("%.6f", cor(rescue_ok$nUMI, spe_cd$sum_umi)), NA, FALSE)
add_qc(
  "nUMI_over_sum_umi_ratio_quantiles",
  paste(sprintf("%.4f", quantile(umi_ratio, c(0, 0.01, 0.5, 0.99, 1))), collapse = ";"),
  NA, FALSE, "min;1%;median;99%;max"
)
add_qc("nUMI_greater_than_sum_umi_n", sum(umi_ratio > 1), sum(umi_ratio > 1) == 0, FALSE)

## 5) Validity of the proportions ----
prop_mat <- as.matrix(rescue_df[, ct_cols])
prop_sum <- rowSums(prop_mat)

add_qc("n_NA_proportions", sum(is.na(prop_mat)), !anyNA(prop_mat), TRUE)
add_qc(
  "proportions_within_0_1", all(prop_mat >= 0 & prop_mat <= 1, na.rm = TRUE),
  all(prop_mat >= 0 & prop_mat <= 1, na.rm = TRUE), TRUE
)
add_qc(
  "max_abs_deviation_rowsum_from_1",
  sprintf("%.3e", max(abs(prop_sum - 1))),
  max(abs(prop_sum - 1)) < prop_sum_tol, TRUE,
  "row sum over all cell types incl. 'Unattributed Expression'"
)

## 6) Two "unattributed" quantities (diagnostic only) ----
# `unattributed_ratio` and the `Unattributed Expression` proportion are
# reported separately by RESCUE and are not identical; both are kept.
unattr_diff <- rescue_df$unattributed_ratio - rescue_df$`Unattributed Expression`
add_qc(
  "unattributed_ratio_minus_Unattributed_Expression_quantiles",
  paste(sprintf("%.4f", quantile(unattr_diff, c(0, 0.01, 0.5, 0.99, 1))), collapse = ";"),
  NA, FALSE, "min;1%;median;99%;max"
)
add_qc(
  "unattributed_ratio_cor_Unattributed_Expression",
  sprintf("%.4f", cor(rescue_df$unattributed_ratio, rescue_df$`Unattributed Expression`)),
  NA, FALSE
)

qc_df <- bind_rows(qc_list)
print(qc_df, n = Inf, width = Inf)

## Per-sample spot coverage ----
qc_sample_df <- full_join(
  colData(spe) |>
    as.data.frame() |>
    count(sample_id, brnum, dx, name = "n_spots_spe"),
  rescue_df |>
    count(sample_id, name = "n_spots_rescue"),
  by = "sample_id"
) |>
  left_join(
    tibble(sample_id = spe$sample_id, matched = spe$key %in% rescue_df$spot_id) |>
      group_by(sample_id) |>
      summarize(n_spots_matched = sum(matched), .groups = "drop"),
    by = "sample_id"
  ) |>
  mutate(
    across(starts_with("n_spots"), \(x) replace_na(x, 0L)),
    frac_spe_matched = n_spots_matched / n_spots_spe
  ) |>
  arrange(frac_spe_matched, sample_id) |>
  as_tibble()

print(qc_sample_df, n = Inf)

write_csv(qc_df, file.path(fld_out, "RESCUE_merge_qc_summary.csv"))
write_csv(qc_sample_df, file.path(fld_out, "RESCUE_merge_qc_per_sample.csv"))

# error prevention: stop on any failed misalignment check
failed_hard <- qc_df |> filter(hard, !pass)
if (nrow(failed_hard) > 0) {
  print(failed_hard)
  stop("RESCUE merge failed hard QC check(s); see RESCUE_merge_qc_summary.csv")
}


# Build slim, key-indexed proportion data frame -----------------------------
# Cell type columns are prefixed with "RESCUE_" and made syntactic
# ("Unattributed Expression" -> "RESCUE_Unattributed").
ct_rename <- set_names(
  ct_cols,
  paste0("RESCUE_", if_else(ct_cols == "Unattributed Expression", "Unattributed", ct_cols))
)

rescue_prop_df <- rescue_df |>
  select(key = spot_id, all_of(ct_cols), nUMI, unattributed_ratio) |>
  rename(
    all_of(ct_rename),
    RESCUE_nUMI = nUMI,
    RESCUE_unattributed_ratio = unattributed_ratio
  )

rescue_vars <- setdiff(colnames(rescue_prop_df), "key")


# Merge into colData(spe) ------------------------------------------------------
# Same convention as PRECAST labels in
# `code/analysis/01_build_spe/20_create_final_spe.R`: match on `key`, keeping
# spe column order. Spots without a RESCUE result (if any) get NA.
m <- match(spe$key, rescue_prop_df$key)
for (.var in rescue_vars) {
  spe[[.var]] <- rescue_prop_df[[.var]][m]
}

# error prevention
stopifnot(all(rescue_vars %in% colnames(colData(spe))))
stopifnot(sum(!is.na(spe$RESCUE_Astro)) == n_spe_in_rescue)
stopifnot(all.equal(
  sum(spe$RESCUE_Excit_01, na.rm = TRUE),
  sum(rescue_df$Excit_01[rescue_df$spot_id %in% spe$key])
))

# Record cell type columns for downstream scripts
metadata(spe)$RESCUE_celltype_cols <- names(ct_rename)

# Long-form slim df with spot metadata for lightweight downstream plotting
rescue_prop_df <- colData(spe) |>
  as.data.frame() |>
  select(key, sample_id, brnum, dx, spd_label, all_of(rescue_vars)) |>
  as_tibble()


# Save -------------------------------------------------------------------------
saveRDS(
  rescue_prop_df,
  file.path(fld_out, "RESCUE_prop_df.rds")
)

saveRDS(
  spe,
  file.path(fld_out, "fnl_spe_kept_spots_only_w_RESCUE.rds")
)

print("Finished merging RESCUE cell type proportions into spe")


# Session info ----
sessioninfo::session_info()
