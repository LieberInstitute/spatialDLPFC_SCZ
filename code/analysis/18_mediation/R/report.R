## tables and base-R figures avoid additional report dependencies.
report_stage <- function(screens, cfg, env) {
  r <- read_tsv(file.path(env$out, "primary", "all_pairs.tsv.gz"))
  assert(identical(unique(r$run_signature), env$signature), "Combined primary table is stale")
  for (id in screens$screen_id) invisible(get_artifact(id, env))
  assert(setequal(unique(r$screen_id), screens$screen_id), "Report screens do not match primary family")
  summary <- do.call(rbind, lapply(screens$screen_id, function(id) {
    z <- r[r$screen_id == id, ]
    data.frame(screen_id = id, observations = z$n_samples[1], donors = z$n_donors[1], tested_genes = nrow(z),
      a_beta = z$a_beta[1], a_p = z$a_p[1], a_q = z$a_q[1], a_gate = z$mediator_gate[1],
      historical_eligible = sum(z$historical_gate), matched_baseline_eligible = sum(z$historical_gate & z$baseline_gate),
      conditional_M_FDR = sum(z$historical_gate & z$baseline_gate & z$mediator_outcome_gate),
      screening_hits = sum(z$screen_hit), hits_with_shrinkage = sum(z$screen_hit & z$coefficient_shrinkage),
      higher_priority_hits = sum(z$higher_priority))
  }))
  write_tsv(summary, file.path(env$out, "report", "screen_summary.tsv"))
  current_status <- function(stage) {
    p <- file.path(env$out, stage, "status.tsv")
    if (!file.exists(p)) return(NULL)
    s <- read_tsv(p)
    if (!identical(unique(s$run_signature), env$signature)) return(NULL)
    s
  }
  sensitivity_state <- current_status("sensitivity")
  overlap_state <- current_status("overlap")
  hits <- r[r$screen_hit, ]; sensitivity_flags <- list()
  for (id in screens$screen_id) {
    z <- hits[hits$screen_id == id, ]
    if (!nrow(z)) next
    files <- list.files(file.path(env$out, "sensitivity"), pattern = paste0("^", id, "_.*[.]tsv[.]gz$"), full.names = TRUE)
    overlap_path <- file.path(env$out, "overlap", paste0(id, "_results.tsv.gz"))
    if (is.null(sensitivity_state)) files <- character()
    if (!is.null(overlap_state) && file.exists(overlap_path)) files <- c(files, overlap_path)
    for (p in files) {
      b <- read_tsv(p)
      assert(identical(unique(b$run_signature), env$signature), paste("Stale sensitivity table:", p))
      ix <- match(z$gene_id, b$gene_id)
      sensitivity_flags[[length(sensitivity_flags) + 1L]] <- data.frame(screen_id = id, gene_id = z$gene_id,
        analysis = b$analysis[1], testable = !is.na(ix), screen_hit = b$screen_hit[ix],
        coefficient_shrinkage = b$coefficient_shrinkage[ix], c_beta = b$c_beta[ix], cprime_beta = b$cprime_beta[ix],
        a_p = b$a_p[ix], b_q = b$b_q[ix], cprime_p = b$cprime_p[ix])
    }
  }
  if (length(sensitivity_flags)) write_tsv(do.call(rbind, sensitivity_flags), file.path(env$out, "report", "hit_sensitivity_comparison.tsv"))
  dir.create(env$plots, recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file.path(env$plots, "primary_diagnostics.pdf"), width = 10, height = 5.5)
  for (id in screens$screen_id) {
    a <- get_artifact(id, env); z <- r[r$screen_id == id, ]; d <- a$samples
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    boxplot(M_logCPM ~ Dx, d, xlab = "Diagnosis", ylab = "Source mediator logCPM", main = id)
    stripchart(M_logCPM ~ Dx, d, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, col = "#33333366")
    mtext(paste(nrow(d), "donor-SpD observations;", length(unique(d$donor)), "donors"), side = 3, cex = 0.7)
    use <- z$historical_gate & z$baseline_gate
    if (any(use)) {
    plot(z$c_beta[use], z$cprime_beta[use], pch = 16, cex = 0.6,
         col = ifelse(z$screen_hit[use], "#C1482F", "#33333355"),
         xlab = "Baseline SCZ - NTC coefficient", ylab = "Mediator-adjusted coefficient", main = "Eligible outcome genes")
    abline(0, 1, lty = 2); abline(h = 0, v = 0, col = "grey80")
    } else { plot.new(); title("No eligible baseline outcome genes") }
  }
  grDevices::dev.off()
  if (!nrow(hits)) {
    unlink(file.path(env$out, "report", "hit_sensitivity_comparison.tsv"))
    unlink(file.path(env$plots, "hit_effects_and_partial_associations.pdf"))
  }
  if (nrow(hits)) {
    grDevices::pdf(file.path(env$plots, "hit_effects_and_partial_associations.pdf"), width = 10, height = 5.5)
    for (id in screens$screen_id) {
      a <- get_artifact(id, env); z <- hits[hits$screen_id == id, ]
      if (!nrow(z)) next
      z <- head(z[order(!z$higher_priority, z$b_q), ], 30)
      d <- a$samples
      donor_cols <- setNames(grDevices::hcl.colors(length(unique(d$donor)), "Dynamic"), unique(d$donor))
      for (i in seq_len(nrow(z))) {
        gene <- z$gene_id[i]
        par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
        values <- c(z$c_beta[i], z$cprime_beta[i]); lo <- c(z$c_lower[i], z$cprime_lower[i]); hi <- c(z$c_upper[i], z$cprime_upper[i])
        plot(values, 1:2, xlim = range(lo, hi), ylim = c(0.5, 2.5), yaxt = "n", pch = 19,
             xlab = "Diagnosis coefficient (95% CI)", ylab = "", main = paste(id, z$gene_name[i]))
        axis(2, at = 1:2, labels = c("Baseline", "+ mediator"), las = 1)
        segments(lo, 1:2, hi, 1:2); abline(v = 0, lty = 2)
        my <- a$baseline$EList$E[gene, ]
        rx <- lm.fit(a$design, d$M)$residuals; ry <- lm.fit(a$design, my)$residuals
        plot(rx, ry, col = donor_cols[d$donor], pch = 19, xlab = "Mediator residual", ylab = "Outcome residual",
             main = "Unweighted partial-association diagnostic")
        abline(lm(ry ~ rx), col = "grey30")
        mtext("Color identifies donor; observations are repeated SpDs", side = 3, cex = 0.65)
      }
    }
    grDevices::dev.off()
  }
  sensitivity_status <- file.path(env$out, "sensitivity", "status.tsv")
  overlap_status <- file.path(env$out, "overlap", "status.tsv")
  sens_text <- if (!is.null(sensitivity_state)) {
    s <- sensitivity_state; paste(sum(s$status == "complete"), "of", nrow(s), "required sensitivity tasks completed.")
  } else "Sensitivity analyses have not run."
  overlap_text <- if (!is.null(overlap_state)) {
    s <- overlap_state; paste("Overlap audit completed.", sum(s$status == "complete"), "screens reaggregated;",
      sum(s$status == "not_required_no_primary_hits"), "did not require reaggregation;",
      sum(s$status == "untestable_or_failed"), "were untestable or failed.")
  } else "Shared-spot audit has not run; overlap robustness is unverified."
  header <- "| Screen | Donors | Genes tested | Mediator Dx p | Hits | Higher priority |"
  rows <- vapply(seq_len(nrow(summary)), function(i) {
    z <- summary[i, ]; sprintf("| %s | %d | %d | %.4g | %d | %d |", z$screen_id, z$donors, z$tested_genes, z$a_p, z$screening_hits, z$higher_priority_hits)
  }, "")
  text <- c("# SCZ PTN/FGF exploratory mediation screening", "", paste("Generated:", Sys.time()),
    paste("Run signature:", env$signature), paste("Model engine:", cfg$engine), "", header, "|---|---:|---:|---:|---:|---:|", rows, "",
    sens_text, overlap_text, "",
    sprintf("Screening rule: historical and matched-baseline diagnosis p < %g, %sconditional mediator-outcome BH FDR < %g, and adjusted diagnosis p >= %g.",
            cfg$dx_p, if (isFALSE(cfg$require_mediator_gate)) "" else sprintf("matched mediator diagnosis p < %g, ", cfg$mediator_dx_p), cfg$mediator_q, cfg$dx_p), "",
    sprintf("Higher priority additionally requires same-direction coefficient shrinkage and BH FDR < %g across all five mediator-outcome testing families. This is not a mediation-discovery FDR guarantee.", cfg$mediator_q), "",
    "All four historical mediator nominations fail gene-wide FDR 0.05. Nominal selection and testing reuse the same cohort. A p-value crossing 0.05 is not a test of coefficient change or proof of mediation.", "",
    "These are cross-sectional mixed-tissue observations with repeated SpDs per donor and potentially shared spots between SPGs. Diagnosis is not randomized; temporal ordering and unmeasured confounding remain unresolved. Results support hypotheses for follow-up, not causal or functional validation.", "",
    "See screen_summary.tsv, hit_sensitivity_comparison.tsv (when hits exist), primary/all_pairs.tsv.gz, historical/reconciliation.tsv, and the sensitivity/overlap status tables for complete evidence.")
  writeLines(text, file.path(env$out, "report", "RESULTS.md"))
}
