#!/bin/env Rscript
suppressMessages({
  library(here)
  library(data.table)
  library(qs2)
  library(coloc)
})

here::i_am('.git/HEAD')

datadir <- here("processed-data", "eQTL")
primary_window <- "w500k"
rule_txt <- "H4 > 0.8 | ((H3+H4) > 0.8 & H4/(H3+H4) > 0.9)"
npoints_cfg <- 200L
p12_plausible_min_cfg <- 5e-6
p12_plausible_max_cfg <- 5e-5
gate_frac_min_cfg <- 0.5

coloc_dir <- file.path(datadir, "coloc", primary_window)
if (!dir.exists(coloc_dir)) stop("coloc directory not found: ", coloc_dir)

f_primary <- file.path(datadir, "coloc", "coloc_abf_results.tsv.gz")
if (!file.exists(f_primary)) {
  stop("Primary coloc table not found: ", f_primary, "\nRun 05_coloc_explore.Rmd first.")
}

assign_cat_from_pp <- function(pp3, pp4) {
  pp34 <- pp3 + pp4
  ratio <- fifelse(!is.na(pp34) & pp34 > 0, pp4 / pp34, NA_real_)
  fcase(
    pp4 > 0.8, "strong_coloc",
    pp34 > 0.8 & ratio >= 0.9, "likely_coloc",
    pp34 > 0.8 & ratio < 0.5, "distinct_causal",
    (pp34 > 0.8) | (ratio > 0.8), "follow_up",
    default = NA_character_
  )
}

primary <- fread(f_primary)
if (!all(c("ds", "fid", "cat", "gene_name", "lead_snp") %chin% names(primary))) {
  stop("Primary coloc table missing required columns: ds, fid, cat, gene_name, lead_snp")
}
primary[, cat_raw := as.character(cat)]
primary[is.na(cat_raw) | cat_raw == "", cat_raw := NA_character_]

target <- unique(
  primary[!is.na(cat_raw), .(ds, fid, cat_raw, gene_name, lead_snp, PP3, PP4, PP34, PP4_over_PP34)],
  by = c("ds", "fid")
)
if (!nrow(target)) stop("No categorized rows found in: ", f_primary)

rows <- list()
for (ctx in sort(unique(target$ds))) {
  f_qs <- file.path(coloc_dir, sprintf("coloc_%s_gene.qs2", ctx))
  if (!file.exists(f_qs)) {
    warning("Skipping ds=", ctx, " (missing file: ", f_qs, ")")
    next
  }
  message(Sys.time(), " | loading ", f_qs)
  l <- qs_read(f_qs)

  sub <- target[ds == ctx]
  for (i in seq_len(nrow(sub))) {
    fid_val <- sub$fid[i]
    cat_val <- sub$cat_raw[i]
    obj <- l[[fid_val]]

    if (is.null(obj) || is.null(obj$summary) || is.null(obj$priors)) {
      rows[[length(rows) + 1L]] <- data.table(
        window = primary_window, ds = ctx, fid = fid_val, cat_raw = cat_val,
        gene_name = sub$gene_name[i], lead_snp = sub$lead_snp[i],
        cat_default = NA_character_, default_cat_match = NA,
        p1 = NA_real_, p2 = NA_real_, p12_default = NA_real_,
        PP3 = sub$PP3[i], PP4 = sub$PP4[i], PP34 = sub$PP34[i], PP4_over_PP34 = sub$PP4_over_PP34[i],
        default_pass = NA, npoints = NA_integer_, n_pass = NA_integer_,
        frac_pass = NA_real_, min_pass_p12 = NA_real_, max_pass_p12 = NA_real_,
        n_plausible = NA_integer_, n_match_cat = NA_integer_, frac_match_cat = NA_real_,
        p12_plausible_min = p12_plausible_min_cfg, p12_plausible_max = p12_plausible_max_cfg,
        rule = rule_txt, error = "missing coloc object/summary/priors"
      )
      next
    }

    sm <- as.list(obj$summary)
    pp3 <- as.numeric(sm[["PP.H3.abf"]])
    pp4 <- as.numeric(sm[["PP.H4.abf"]])
    pp34 <- pp3 + pp4
    ratio <- if (is.na(pp34) || pp34 <= 0) NA_real_ else pp4 / pp34
    cat_default_val <- assign_cat_from_pp(pp3, pp4)
    default_cat_match_val <- !is.na(cat_default_val) && cat_default_val == cat_val
    default_pass_val <- default_cat_match_val

    sens <- tryCatch(
      coloc::sensitivity(obj, rule = rule_txt, doplot = FALSE, npoints = npoints_cfg),
      error = function(e) e
    )

    if (inherits(sens, "error")) {
      rows[[length(rows) + 1L]] <- data.table(
        window = primary_window, ds = ctx, fid = fid_val, cat_raw = cat_val,
        gene_name = sub$gene_name[i], lead_snp = sub$lead_snp[i],
        cat_default = cat_default_val, default_cat_match = default_cat_match_val,
        p1 = as.numeric(obj$priors[["p1"]]),
        p2 = as.numeric(obj$priors[["p2"]]),
        p12_default = as.numeric(obj$priors[["p12"]]),
        PP3 = pp3, PP4 = pp4, PP34 = pp34, PP4_over_PP34 = ratio,
        default_pass = default_pass_val, npoints = NA_integer_, n_pass = NA_integer_,
        frac_pass = NA_real_, min_pass_p12 = NA_real_, max_pass_p12 = NA_real_,
        n_plausible = NA_integer_, n_match_cat = NA_integer_, frac_match_cat = NA_real_,
        p12_plausible_min = p12_plausible_min_cfg, p12_plausible_max = p12_plausible_max_cfg,
        rule = rule_txt, error = conditionMessage(sens)
      )
      next
    }

    sens_dt <- as.data.table(sens)
    sens_dt[, cat_at_p12 := assign_cat_from_pp(`PP.H3.abf`, `PP.H4.abf`)]
    n_pass_val <- sens_dt[, sum(pass, na.rm = TRUE)]
    n_tot <- nrow(sens_dt)
    pass_p12 <- sens_dt[pass == TRUE, p12]
    min_pass <- if (length(pass_p12)) min(pass_p12, na.rm = TRUE) else NA_real_
    max_pass <- if (length(pass_p12)) max(pass_p12, na.rm = TRUE) else NA_real_
    plausible <- sens_dt[p12 >= p12_plausible_min_cfg & p12 <= p12_plausible_max_cfg]
    n_plausible_val <- nrow(plausible)
    n_match_cat_val <- if (n_plausible_val > 0L) plausible[, sum(cat_at_p12 == cat_val, na.rm = TRUE)] else 0L
    frac_match_cat_val <- if (n_plausible_val > 0L) n_match_cat_val / n_plausible_val else NA_real_

    rows[[length(rows) + 1L]] <- data.table(
      window = primary_window, ds = ctx, fid = fid_val, cat_raw = cat_val,
      gene_name = sub$gene_name[i], lead_snp = sub$lead_snp[i],
      cat_default = cat_default_val, default_cat_match = default_cat_match_val,
      p1 = as.numeric(obj$priors[["p1"]]),
      p2 = as.numeric(obj$priors[["p2"]]),
      p12_default = as.numeric(obj$priors[["p12"]]),
      PP3 = pp3, PP4 = pp4, PP34 = pp34, PP4_over_PP34 = ratio,
      default_pass = default_pass_val, npoints = n_tot, n_pass = n_pass_val,
      frac_pass = if (n_tot > 0) n_pass_val / n_tot else NA_real_,
      min_pass_p12 = min_pass, max_pass_p12 = max_pass,
      n_plausible = n_plausible_val, n_match_cat = n_match_cat_val, frac_match_cat = frac_match_cat_val,
      p12_plausible_min = p12_plausible_min_cfg, p12_plausible_max = p12_plausible_max_cfg,
      rule = rule_txt, error = NA_character_
    )
  }
}

if (!length(rows)) stop("No prior sensitivity rows were generated.")
res <- rbindlist(rows, fill = TRUE)
res[, cat := cat_raw]

res[, gate_pass := fcase(
  is.na(cat_raw), NA,
  !is.na(error), FALSE,
  is.na(default_cat_match) | default_cat_match == FALSE, FALSE,
  is.na(n_plausible) | n_plausible <= 0L, FALSE,
  is.na(frac_match_cat), FALSE,
  frac_match_cat < gate_frac_min_cfg, FALSE,
  default = TRUE
)]

res[, gate_reason := fcase(
  is.na(cat_raw), "not_evaluated",
  !is.na(error), "fail_sensitivity_error",
  is.na(default_cat_match) | default_cat_match == FALSE, "fail_default_rule",
  is.na(n_plausible) | n_plausible <= 0L, "fail_no_plausible_points",
  is.na(frac_match_cat) | frac_match_cat < gate_frac_min_cfg, "fail_cat_instability",
  default = "pass"
)]
res[, cat_gated := fifelse(gate_pass == TRUE, cat_raw, NA_character_)]
res[, gate_frac_min := gate_frac_min_cfg]

summary_dt <- res[, .(
  n_tested = .N,
  n_with_sensitivity = sum(!is.na(npoints)),
  n_default_pass = sum(default_pass == TRUE, na.rm = TRUE),
  n_default_fail = sum(default_pass == FALSE, na.rm = TRUE),
  n_any_pass_grid = sum(!is.na(min_pass_p12)),
  n_no_pass_grid = sum(is.na(min_pass_p12)),
  n_error = sum(!is.na(error)),
  n_gate_pass = sum(gate_pass == TRUE, na.rm = TRUE),
  n_gate_fail = sum(gate_pass == FALSE, na.rm = TRUE),
  default_pass_rate = mean(default_pass == TRUE, na.rm = TRUE),
  any_pass_grid_rate = mean(!is.na(min_pass_p12)),
  gate_pass_rate = mean(gate_pass == TRUE, na.rm = TRUE)
), by = .(window, ds, cat_raw)][order(ds, cat_raw)]
summary_dt[, cat := cat_raw]

gated_cols <- c(
  "window", "ds", "fid", "cat_raw", "cat_default", "default_cat_match",
  "p1", "p2", "p12_default", "p12_plausible_min", "p12_plausible_max",
  "npoints", "n_pass", "frac_pass", "min_pass_p12", "max_pass_p12",
  "n_plausible", "n_match_cat", "frac_match_cat",
  "gate_frac_min", "gate_pass", "gate_reason", "cat_gated",
  "rule", "error"
)
gated_core <- res[, ..gated_cols]
primary_base <- copy(primary)
drop_overlap <- intersect(c("window", "p1", "p2", "cat_raw"), names(primary_base))
if (length(drop_overlap) > 0L) primary_base[, (drop_overlap) := NULL]
gated <- merge(primary_base, gated_core, by = c("ds", "fid"), all.x = TRUE, sort = FALSE)
gated[is.na(cat_raw), cat_raw := as.character(cat)]
gated[is.na(cat_raw) | cat_raw == "", cat_raw := NA_character_]
gated[is.na(cat_raw), `:=`(
  gate_reason = "not_evaluated",
  gate_pass = NA,
  cat_gated = NA_character_
)]
gated[!is.na(cat_raw) & (is.na(gate_reason) | gate_reason == ""), `:=`(
  gate_reason = "fail_sensitivity_error",
  gate_pass = FALSE,
  cat_gated = NA_character_,
  error = fifelse(is.na(error), "missing sensitivity row", error)
)]
setcolorder(gated, c(
  "window", "ds", "fid", "gene_name", "cat", "cat_raw", "cat_default", "cat_gated",
  "gate_pass", "gate_reason", "default_cat_match", "frac_match_cat", "n_plausible",
  setdiff(names(gated), c(
    "window", "ds", "fid", "gene_name", "cat", "cat_raw", "cat_default", "cat_gated",
    "gate_pass", "gate_reason", "default_cat_match", "frac_match_cat", "n_plausible"
  ))
))
req_gate_cols <- c("cat_raw", "cat_gated", "gate_pass", "gate_reason", "frac_match_cat", "n_plausible")
if (!all(req_gate_cols %chin% names(gated))) {
  stop("gated table missing required columns: ", paste(setdiff(req_gate_cols, names(gated)), collapse = ", "))
}

gated_eval <- gated[!is.na(cat_raw)]
gated_eval[, window := primary_window]
gated_lean_cols <- c(
  "window", "ds", "fid", "gene_name", "lead_snp", "rsid", "nsnps",
  "PP3", "PP4", "PP34", "PP4_over_PP34", "cat_raw", "cat_gated",
  "gate_pass", "gate_reason"
)
if (!all(gated_lean_cols %chin% names(gated_eval))) {
  stop("gated evaluated table missing required columns: ", paste(setdiff(gated_lean_cols, names(gated_eval)), collapse = ", "))
}
gated_eval_lean <- gated_eval[, ..gated_lean_cols]
if (gated_eval_lean[is.na(window), .N] > 0L) stop("gated evaluated table has missing window values")
if (gated_eval_lean[is.na(cat_raw), .N] > 0L) stop("gated evaluated table has missing cat_raw values")
if (gated_eval_lean[is.na(gate_reason) | gate_reason == "", .N] > 0L) stop("gated evaluated table has missing gate_reason values")
if (gated_eval_lean[gate_reason == "not_evaluated", .N] > 0L) stop("gated evaluated table should not contain not_evaluated rows")
if (gated_eval_lean[gate_pass == TRUE & is.na(cat_gated), .N] > 0L) stop("gate_pass/cat_gated mismatch for pass rows")
if (gated_eval_lean[gate_pass == FALSE & !is.na(cat_gated), .N] > 0L) stop("gate_pass/cat_gated mismatch for fail rows")

fail_eval <- copy(gated_eval)
if (fail_eval[is.na(gate_reason), .N] > 0L) stop("categorized rows with missing gate_reason detected")
if (fail_eval[gate_reason == "pass" & gate_pass != TRUE, .N] > 0L) stop("gate_reason/gate_pass mismatch for pass rows")
if (fail_eval[gate_pass == TRUE & gate_reason != "pass", .N] > 0L) stop("gate_reason/gate_pass mismatch for failing rows")

fail_summary <- fail_eval[, .(n_loci = .N), by = .(window, ds, cat_raw, gate_reason)][
  order(ds, cat_raw, gate_reason)
]
fail_summary[, pct_within_ds_cat := n_loci / sum(n_loci), by = .(window, ds, cat_raw)]

fail_summary_global <- fail_eval[, .(n_loci = .N), by = .(window, cat_raw, gate_reason)][
  order(cat_raw, gate_reason)
]
fail_summary_global[, pct_global := n_loci / sum(n_loci), by = .(window)]
fail_table <- fail_summary[window == primary_window][order(ds, cat_raw, gate_reason)]

chk_local <- fail_summary[, .(n_sum = sum(n_loci)), by = .(window, ds, cat_raw)]
chk_base <- fail_eval[, .(n_eval = .N), by = .(window, ds, cat_raw)]
chk <- merge(chk_local, chk_base, by = c("window", "ds", "cat_raw"), all = TRUE)
if (chk[is.na(n_sum) | is.na(n_eval) | n_sum != n_eval, .N] > 0L) {
  stop("fail_summary per context/category does not match evaluated row counts")
}
if (sum(fail_summary$n_loci) != nrow(fail_eval)) {
  stop("global fail_summary total does not match evaluated row count")
}

out_dir <- file.path(datadir, "coloc", "prior_sensitivity")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

f_per_gene <- file.path(out_dir, sprintf("%s_prior_sensitivity_per_gene.tsv.gz", primary_window))
f_summary <- file.path(out_dir, sprintf("%s_prior_sensitivity_summary.tsv.gz", primary_window))
f_fail_summary <- file.path(out_dir, sprintf("%s_prior_sensitivity_fail_summary.tsv.gz", primary_window))
f_fail_summary_global <- file.path(out_dir, sprintf("%s_prior_sensitivity_fail_summary_global.tsv.gz", primary_window))
f_fail_table <- file.path(out_dir, sprintf("%s_prior_sensitivity_fail_table.tsv.gz", primary_window))
f_gated <- file.path(datadir, "coloc", "coloc_abf_results_gated.tsv.gz")

fwrite(res, f_per_gene, sep = "\t", na = "NA", quote = FALSE)
fwrite(summary_dt, f_summary, sep = "\t", na = "NA", quote = FALSE)
fwrite(fail_summary, f_fail_summary, sep = "\t", na = "NA", quote = FALSE)
fwrite(fail_summary_global, f_fail_summary_global, sep = "\t", na = "NA", quote = FALSE)
fwrite(fail_table, f_fail_table, sep = "\t", na = "NA", quote = FALSE)
fwrite(gated_eval_lean, f_gated, sep = "\t", na = "NA", quote = FALSE)

message("Saved per-gene prior sensitivity: ", f_per_gene)
message("Saved summary prior sensitivity: ", f_summary)
message("Saved fail summary: ", f_fail_summary)
message("Saved global fail summary: ", f_fail_summary_global)
message("Saved fail table: ", f_fail_table)
message("Saved gated coloc results (evaluated rows, lean schema): ", f_gated)
