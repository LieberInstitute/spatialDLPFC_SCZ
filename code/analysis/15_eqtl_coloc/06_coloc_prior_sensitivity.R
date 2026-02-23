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
npoints <- 200L

coloc_dir <- file.path(datadir, "coloc", primary_window)
if (!dir.exists(coloc_dir)) stop("coloc directory not found: ", coloc_dir)

f_primary <- file.path(datadir, "coloc", "coloc_abf_results.tsv.gz")
if (!file.exists(f_primary)) {
  stop("Primary coloc table not found: ", f_primary, "\nRun 05_coloc_explore.Rmd first.")
}

target <- fread(f_primary)
target <- target[
  cat %chin% c("strong_coloc", "likely_coloc"),
  .(ds, fid, cat, gene_name, lead_snp, PP3, PP4, PP34, PP4_over_PP34)
]
if (!nrow(target)) stop("No strong_coloc or likely_coloc rows found in: ", f_primary)

eval_rule_default <- function(pp3, pp4) {
  pp34 <- pp3 + pp4
  ratio <- if (is.na(pp34) || pp34 <= 0) NA_real_ else pp4 / pp34
  as.logical((pp4 > 0.8) || (pp34 > 0.8 && ratio > 0.9))
}

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
    fid <- sub$fid[i]
    cat_i <- sub$cat[i]
    obj <- l[[fid]]

    if (is.null(obj) || is.null(obj$summary) || is.null(obj$priors)) {
      rows[[length(rows) + 1L]] <- data.table(
        window = primary_window, ds = ctx, fid = fid, cat = cat_i,
        gene_name = sub$gene_name[i], lead_snp = sub$lead_snp[i],
        p1 = NA_real_, p2 = NA_real_, p12_default = NA_real_,
        PP3 = sub$PP3[i], PP4 = sub$PP4[i], PP34 = sub$PP34[i], PP4_over_PP34 = sub$PP4_over_PP34[i],
        default_pass = NA, npoints = NA_integer_, n_pass = NA_integer_,
        frac_pass = NA_real_, min_pass_p12 = NA_real_, max_pass_p12 = NA_real_,
        rule = rule_txt, error = "missing coloc object/summary/priors"
      )
      next
    }

    sm <- as.list(obj$summary)
    pp3 <- as.numeric(sm[["PP.H3.abf"]])
    pp4 <- as.numeric(sm[["PP.H4.abf"]])
    pp34 <- pp3 + pp4
    ratio <- if (is.na(pp34) || pp34 <= 0) NA_real_ else pp4 / pp34
    default_pass <- eval_rule_default(pp3, pp4)

    sens <- tryCatch(
      coloc::sensitivity(obj, rule = rule_txt, doplot = FALSE, npoints = npoints),
      error = function(e) e
    )

    if (inherits(sens, "error")) {
      rows[[length(rows) + 1L]] <- data.table(
        window = primary_window, ds = ctx, fid = fid, cat = cat_i,
        gene_name = sub$gene_name[i], lead_snp = sub$lead_snp[i],
        p1 = as.numeric(obj$priors[["p1"]]),
        p2 = as.numeric(obj$priors[["p2"]]),
        p12_default = as.numeric(obj$priors[["p12"]]),
        PP3 = pp3, PP4 = pp4, PP34 = pp34, PP4_over_PP34 = ratio,
        default_pass = default_pass, npoints = NA_integer_, n_pass = NA_integer_,
        frac_pass = NA_real_, min_pass_p12 = NA_real_, max_pass_p12 = NA_real_,
        rule = rule_txt, error = conditionMessage(sens)
      )
      next
    }

    sens_dt <- as.data.table(sens)
    n_pass <- sens_dt[, sum(pass, na.rm = TRUE)]
    n_tot <- nrow(sens_dt)
    pass_p12 <- sens_dt[pass == TRUE, p12]
    min_pass <- if (length(pass_p12)) min(pass_p12, na.rm = TRUE) else NA_real_
    max_pass <- if (length(pass_p12)) max(pass_p12, na.rm = TRUE) else NA_real_

    rows[[length(rows) + 1L]] <- data.table(
      window = primary_window, ds = ctx, fid = fid, cat = cat_i,
      gene_name = sub$gene_name[i], lead_snp = sub$lead_snp[i],
      p1 = as.numeric(obj$priors[["p1"]]),
      p2 = as.numeric(obj$priors[["p2"]]),
      p12_default = as.numeric(obj$priors[["p12"]]),
      PP3 = pp3, PP4 = pp4, PP34 = pp34, PP4_over_PP34 = ratio,
      default_pass = default_pass, npoints = n_tot, n_pass = n_pass,
      frac_pass = if (n_tot > 0) n_pass / n_tot else NA_real_,
      min_pass_p12 = min_pass, max_pass_p12 = max_pass,
      rule = rule_txt, error = NA_character_
    )
  }
}

if (!length(rows)) stop("No prior sensitivity rows were generated.")
res <- rbindlist(rows, fill = TRUE)

summary_dt <- res[, .(
  n_tested = .N,
  n_with_sensitivity = sum(!is.na(npoints)),
  n_default_pass = sum(default_pass == TRUE, na.rm = TRUE),
  n_default_fail = sum(default_pass == FALSE, na.rm = TRUE),
  n_any_pass_grid = sum(!is.na(min_pass_p12)),
  n_no_pass_grid = sum(is.na(min_pass_p12)),
  n_error = sum(!is.na(error)),
  default_pass_rate = mean(default_pass == TRUE, na.rm = TRUE),
  any_pass_grid_rate = mean(!is.na(min_pass_p12))
), by = .(window, ds, cat)][order(ds, cat)]

out_dir <- file.path(datadir, "coloc", "prior_sensitivity")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

f_per_gene <- file.path(out_dir, sprintf("%s_prior_sensitivity_per_gene.tsv.gz", primary_window))
f_summary <- file.path(out_dir, sprintf("%s_prior_sensitivity_summary.tsv.gz", primary_window))

fwrite(res, f_per_gene, sep = "\t")
fwrite(summary_dt, f_summary, sep = "\t")

message("Saved per-gene prior sensitivity: ", f_per_gene)
message("Saved summary prior sensitivity: ", f_summary)
