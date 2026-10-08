## audit, primary screens, and deterministic robustness analyses.
audit_stage <- function(inputs, screens, cfg, env) {
  inventory <- list(); matched <- list()
  for (context in names(inputs$objects)) {
    obj <- inputs$objects[[context]]; h <- inputs$historical[[context]]
    exp <- cfg$expected[cfg$expected$context == context, ]
    inventory[[context]] <- data.frame(context = context, genes = nrow(obj$counts),
      samples = ncol(obj$counts), donors = length(unique(obj$samples$donor)),
      historical_genes = nrow(h), historical_nominal = sum(h$p_value_scz < cfg$dx_p),
      historical_fdr05 = sum(h$fdr_scz < 0.05))
    assert(nrow(h) == exp$genes && sum(h$p_value_scz < cfg$dx_p) == exp$nominal,
           paste("Historical audit changed for", context, "- review updated inputs explicitly"))
    write_tsv(obj$samples, file.path(env$out, "audit", paste0(context, "_available_samples.tsv")))
  }
  for (i in seq_len(nrow(screens))) {
    s <- screens[i, ]; src <- inputs$objects[[s$source]]; tgt <- inputs$objects[[s$target]]
    m <- match_samples(src, tgt, cfg); d <- m$samples
    sdge <- prepare_dge(src, d$key, m$design); tdge <- prepare_dge(tgt, d$key, m$design)
    assert(s$mediator_id %in% rownames(sdge), paste(s$screen_id, "mediator unavailable"))
    assert(nrow(d) == s$expected_samples && length(unique(d$donor)) == s$expected_donors,
           paste("Matched sample counts changed for", s$screen_id, "- review configuration"))
    enhanced <- make_design(d, cfg, c("slide_id", "rin"))
    matched[[i]] <- data.frame(screen_id = s$screen_id, observations = nrow(d), donors = length(unique(d$donor)),
      NTC = length(unique(d$donor[d$Dx == "NTC"])), SCZ = length(unique(d$donor[d$Dx == "SCZ"])),
      source_genes = nrow(sdge), target_genes = nrow(tdge), baseline_rank = qr(m$design)$rank,
      enhanced_rank = qr(enhanced)$rank)
    write_tsv(d, file.path(env$out, "audit", paste0(s$screen_id, "_samples.tsv")))
    write_tsv(m$exclusions, file.path(env$out, "audit", paste0(s$screen_id, "_exclusions.tsv")))
    write_tsv(data.frame(gene_id = rownames(tdge)), file.path(env$out, "audit", paste0(s$screen_id, "_target_genes.tsv")))
  }
  write_tsv(do.call(rbind, inventory), file.path(env$out, "audit", "inventory.tsv"))
  write_tsv(do.call(rbind, matched), file.path(env$out, "audit", "matched_designs.tsv"))
}
historical_stage <- function(inputs, cfg, env, workers = 1L) {
  checks <- parallel::mclapply(names(inputs$objects), function(context) {
    obj <- inputs$objects[[context]]
    use <- if (context == "vasc") rep(TRUE, nrow(obj$samples)) else obj$samples$SpD %in% cfg$spds
    d <- obj$samples[use, ]; E <- obj$logcounts[, use, drop = FALSE]
    des <- make_design(d, cfg)
    cd <- make_design(d, cfg, include_spd = FALSE)
    f <- fit_logcounts(E, des, d, cd)
    r <- coef_table(f$fit)
    h <- inputs$historical[[context]]; ix <- match(r$gene_id, h$gene_id)
    r$historical_beta <- h$logFC_scz[ix]; r$historical_p <- h$p_value_scz[ix]
    write_tsv(r, file.path(env$out, "historical", paste0(context, "_reproduction.tsv.gz")))
    delta <- max(abs(r$beta - r$historical_beta), na.rm = TRUE)
    pdiff <- max(abs(r$p - r$historical_p), na.rm = TRUE)
    data.frame(context = context, max_beta_difference = delta,
      max_p_difference = pdiff, reproduced_nominal = sum(r$p < 0.05), historical_nominal = sum(h$p_value_scz < 0.05),
      close_numerical_match = delta < 1e-5 && pdiff < 1e-5)
  }, mc.cores = min(workers, length(inputs$objects)), mc.preschedule = FALSE)
  assert(all(vapply(checks, is.data.frame, logical(1))), "Historical reproduction worker failed")
  write_tsv(do.call(rbind, checks), file.path(env$out, "historical", "reconciliation.tsv"))
}
primary_stage <- function(inputs, screens, all_screens, cfg, env, workers) {
  jobs <- parallel::mclapply(seq_len(nrow(screens)), function(i) {
    s <- screens[i, ]; message("Primary screen: ", s$screen_id)
    tryCatch({
      a <- run_screen(s, inputs, cfg, env)
      path <- file.path(env$out, "primary", paste0(s$screen_id, ".rds"))
      save_atomic(a, path)
      write_tsv(a$samples, file.path(env$out, "primary", paste0(s$screen_id, "_samples.tsv")))
      write_tsv(a$results, file.path(env$out, "primary", paste0(s$screen_id, "_results.tsv.gz")))
      data.frame(screen_id = s$screen_id, status = "complete", message = "", artifact = path,
                 baseline_key = attr(a$baseline, "cache_key"), joint_key = attr(a$joint, "cache_key"))
    }, error = function(e) data.frame(screen_id = s$screen_id, status = "failed", message = conditionMessage(e),
                                      artifact = "", baseline_key = "", joint_key = ""))
  }, mc.cores = min(workers, nrow(screens)), mc.preschedule = FALSE)
  assert(all(vapply(jobs, is.data.frame, logical(1))), "Worker process failure; inspect scheduler log")
  status <- do.call(rbind, jobs)
  write_tsv(status, file.path(env$out, "primary", "status.tsv"))
  assert(all(status$status == "complete"), "One or more screens failed; see primary/status.tsv")
  r <- do.call(rbind, lapply(status$artifact, function(p) readRDS(p)$results))
  complete_family <- setequal(screens$screen_id, all_screens$screen_id)
  if (complete_family) r$b_q_global <- p.adjust(r$b_p, method = "BH")
  r$pooled_family_complete <- complete_family
  r <- classify_pairs(r, cfg)
  r <- r[order(!r$higher_priority, !r$screen_hit, r$b_q_global, r$b_q, -r$absolute_shrinkage), ]
  write_tsv(r, file.path(env$out, "primary", "all_pairs.tsv.gz"))
  write_tsv(r[r$screen_hit, ], file.path(env$out, "primary", "screening_hits.tsv"))
  meds <- unique(r[, c("screen_id", "a_beta", "a_p", "a_q", "a_direction_agrees", "mediator_gate")])
  meds$a_q_five_candidates <- if (complete_family) p.adjust(meds$a_p, "BH") else NA_real_
  write_tsv(meds, file.path(env$out, "primary", "mediator_gates.tsv"))
}
get_artifact <- function(id, env) {
  a <- readRDS(file.path(env$out, "primary", paste0(id, ".rds")))
  assert(identical(a$signature, env$signature), "Primary artifacts are stale; rerun primary screening")
  a
}
sensitivity_result <- function(a, baseline, joint, source_fit, cfg, inputs, label) {
  r <- assemble_results(a$screen, inputs$historical[[a$screen$target]], inputs$objects[[a$screen$target]]$genes,
                         source_fit, baseline, joint, cfg)
  r$analysis <- label; r$run_signature <- a$signature; r$n_samples <- nrow(a$samples); r$n_donors <- length(unique(a$samples$donor))
  r
}
influence_diagnostics <- function(a, ids, cfg, fixed, env) {
  d <- a$samples; E <- a$baseline$EList
  ix <- match(ids, rownames(E$E)); records <- list()
  for (donor in unique(d$donor)) {
    keep <- d$donor != donor; ds <- d[keep, ]
    z <- tryCatch({
      b <- limma::lmFit(E$E[ix, keep, drop = FALSE], make_design(ds, cfg), block = ds$donor,
                       correlation = fixed$rho, weights = E$weights[ix, keep, drop = FALSE])
      j <- limma::lmFit(E$E[ix, keep, drop = FALSE], make_design(ds, cfg, "M"), block = ds$donor,
                       correlation = fixed$rho, weights = E$weights[ix, keep, drop = FALSE])
      data.frame(gene_id = ids, excluded_donor = donor, c_beta = b$coefficients[, "DxSCZ"] - b$coefficients[, "DxNTC"],
        cprime_beta = j$coefficients[, "DxSCZ"] - j$coefficients[, "DxNTC"], b_beta = j$coefficients[, "M"], status = "complete")
    }, error = function(e) data.frame(gene_id = ids, excluded_donor = donor, c_beta = NA_real_, cprime_beta = NA_real_, b_beta = NA_real_, status = conditionMessage(e)))
    records[[donor]] <- z
  }
  write_tsv(do.call(rbind, records), file.path(env$out, "diagnostics", paste0(a$screen$screen_id, "_leave_donor_out.tsv")))
  d$M_between <- ave(d$M, d$donor, FUN = mean)
  d$M_within <- d$M - d$M_between
  des <- make_design(d, cfg, c("M_between", "M_within"))
  f <- limma::lmFit(E$E, des, block = d$donor, correlation = fixed$rho, weights = E$weights)
  between <- prefix_coef(coef_table(f, "M_between"), "between")
  within <- prefix_coef(coef_table(f, "M_within"), "within")
  r <- merge(between, within, by = "gene_id")
  r$primary_hit <- r$gene_id %in% ids
  write_tsv(r, file.path(env$out, "diagnostics", paste0(a$screen$screen_id, "_between_within.tsv.gz")))
}
sensitivity_stage <- function(inputs, screens, cfg, env, workers = 1L) {
  primary <- read_tsv(file.path(env$out, "primary", "all_pairs.tsv.gz"))
  status <- parallel::mclapply(seq_len(nrow(screens)), function(i) {
    status <- list()
    s <- screens[i, ]; a <- get_artifact(s$screen_id, env); d <- a$samples
    message("Sensitivity analyses: ", s$screen_id)
    run <- function(label, fun) {
      tryCatch({
        r <- fun()
        write_tsv(r, file.path(env$out, "sensitivity", paste0(s$screen_id, "_", label, ".tsv.gz")))
        data.frame(screen_id = s$screen_id, analysis = label, status = "complete", message = "")
      }, error = function(e) data.frame(screen_id = s$screen_id, analysis = label, status = "failed", message = conditionMessage(e)))
    }
    status[[length(status) + 1L]] <- run("slide_rin", function() {
      des <- make_design(d, cfg, c("slide_id", "rin"))
      sf <- cached_voom(a$sdge, des, d, env)
      bf <- cached_voom(a$tdge, des, d, env)
      jf <- cached_voom(a$tdge, make_design(d, cfg, c("slide_id", "rin", "M")), d, env)
      sensitivity_result(a, bf, jf, sf, cfg, inputs, "slide_rin")
    })
    status[[length(status) + 1L]] <- run("matched_logcounts", function() {
      ts <- inputs$objects[[s$target]]$logcounts[rownames(a$tdge), d$key, drop = FALSE]
      ss <- inputs$objects[[s$source]]$logcounts[rownames(a$sdge), d$key, drop = FALSE]
      dl <- d; dl$M <- as.numeric(scale(ss[s$mediator_id, ]))
      sf <- fit_logcounts(ss, a$design, dl)$fit
      bf <- fit_logcounts(ts, a$design, dl)$fit
      jf <- fit_logcounts(ts, make_design(dl, cfg, "M"), dl)$fit
      sensitivity_result(a, bf, jf, sf, cfg, inputs, "matched_logcounts")
    })
    fixed <- NULL
    status[[length(status) + 1L]] <- run("fixed_weights", function() {
      fixed <<- fixed_comparison(a, cfg)
      sensitivity_result(a, fixed$baseline, fixed$joint, a$source_fit, cfg, inputs, "fixed_weights")
    })
    if (s$source == "neun" || s$target == "neun") {
      status[[length(status) + 1L]] <- run("without_spd07", function() {
        b <- run_screen(s, inputs, cfg, env, setdiff(cfg$spds, "spd07"),
                        list(source = rownames(a$sdge), target = rownames(a$tdge)))
        b$results$analysis <- "without_spd07"; b$results
      })
    }
    ids <- primary$gene_id[primary$screen_id == s$screen_id & primary$screen_hit]
    if (length(ids)) {
      status[[length(status) + 1L]] <- tryCatch({
        assert(!is.null(fixed), "Fixed-weight fit failed")
        influence_diagnostics(a, ids, cfg, fixed, env)
        data.frame(screen_id = s$screen_id, analysis = "donor_influence_between_within", status = "complete", message = "")
      }, error = function(e) data.frame(screen_id = s$screen_id, analysis = "donor_influence_between_within", status = "failed", message = conditionMessage(e)))
    }
    status
  }, mc.cores = min(workers, nrow(screens)), mc.preschedule = FALSE)
  assert(!any(vapply(status, inherits, logical(1), "try-error")), "Sensitivity worker failed")
  status <- unlist(status, recursive = FALSE)
  fgf_ids <- c("fgf1_neuropil_vasc", "fgf2_neuropil_vasc")
  if (all(fgf_ids %in% screens$screen_id) && any(primary$screen_hit & primary$screen_id %in% fgf_ids)) {
    status[[length(status) + 1L]] <- tryCatch({
      a <- get_artifact(fgf_ids[1], env); b <- get_artifact(fgf_ids[2], env)
      assert(identical(a$samples$key, b$samples$key), "FGF samples differ")
      d <- a$samples; d$M2 <- b$samples$M
      des <- make_design(d, cfg, c("M", "M2"))
      f <- cached_voom(a$tdge, des, d, env)
      r <- merge(prefix_coef(coef_table(f, "M"), "FGF1"), prefix_coef(coef_table(f, "M2"), "FGF2"), by = "gene_id")
      r$mediator_correlation <- cor(d$M, d$M2); r$design_condition_number <- kappa(des)
      write_tsv(r, file.path(env$out, "diagnostics", "neuropil_FGF1_FGF2_joint.tsv.gz"))
      data.frame(screen_id = "neuropil_FGF1_FGF2", analysis = "joint_mediators", status = "complete", message = "")
    }, error = function(e) data.frame(screen_id = "neuropil_FGF1_FGF2", analysis = "joint_mediators", status = "failed", message = conditionMessage(e)))
  }
  status <- do.call(rbind, status)
  status$run_signature <- env$signature
  write_tsv(status, file.path(env$out, "sensitivity", "status.tsv"))
  assert(all(status$status == "complete"), "Sensitivity failure; inspect sensitivity/status.tsv")
}
