## quantify shared spots and rebuild only screens with primary hits.
overlap_stage <- function(inputs, screens, cfg, env) {
  primary <- read_tsv(file.path(env$out, "primary", "all_pairs.tsv.gz"))
  for (id in screens$screen_id) invisible(get_artifact(id, env))
  path <- file.path(env$root, cfg$spots)
  assert(file.exists(path), "Raw spot object unavailable")
  message("Loading raw spot object for overlap audit")
  raw <- readRDS(path)
  cd <- as.data.frame(SummarizedExperiment::colData(raw))
  assert(all(c("sample_id", "fnl_spd", "vasc_pos", "neun_pos", "neuropil_pos") %in% names(cd)), "Raw spot metadata schema mismatch")
  maps <- unique(do.call(rbind, lapply(inputs$objects, function(o) o$samples[, c("donor", "section", "Dx")])))
  unique_keys(maps$section, "Section-to-donor mapping")
  donor <- maps$donor[match(cd$sample_id, maps$section)]
  assert(!anyNA(donor), "Raw section has no donor mapping")
  spd <- canonical_spd(cd$fnl_spd)
  keys <- paste(donor, spd, sep = "|")
  gm <- spd %in% cfg$spds
  keys_gm <- sort(unique(keys[gm])); group <- match(keys, keys_gm)
  labels <- lapply(c(vasc = "vasc_pos", neun = "neun_pos", neuropil = "neuropil_pos"), function(n) {
    v <- cd[[n]]; assert(!anyNA(v), paste("Missing SPG labels:", n)); as.logical(v)
  })
  pair_names <- unique(vapply(seq_len(nrow(screens)), function(i) paste(sort(c(screens$source[i], screens$target[i])), collapse = "|"), ""))
  audit <- list()
  for (pair in pair_names) {
    contexts <- strsplit(pair, "|", fixed = TRUE)[[1]]
    a <- labels[[contexts[1]]]; b <- labels[[contexts[2]]]
    n_a <- tabulate(group[gm & a], nbins = length(keys_gm))
    n_b <- tabulate(group[gm & b], nbins = length(keys_gm))
    shared <- tabulate(group[gm & a & b], nbins = length(keys_gm))
    donor_key <- sub("[|].*$", "", keys_gm)
    audit[[pair]] <- data.frame(pair = pair, key = keys_gm, donor = donor_key,
      SpD = sub("^.*[|]", "", keys_gm), Dx = maps$Dx[match(donor_key, maps$donor)],
      context_a = contexts[1], context_b = contexts[2], n_a = n_a, n_b = n_b, shared = shared,
      fraction_a = ifelse(n_a > 0, shared / n_a, NA_real_), fraction_b = ifelse(n_b > 0, shared / n_b, NA_real_))
  }
  audit <- do.call(rbind, audit)
  write_tsv(audit, file.path(env$out, "overlap", "spot_overlap.tsv"))
  write_tsv(aggregate(cbind(n_a, n_b, shared) ~ pair + Dx, audit, sum), file.path(env$out, "overlap", "overlap_by_diagnosis.tsv"))
  gene_ids <- canonical_gene(SummarizedExperiment::rowData(raw)$gene_id)
  unique_keys(gene_ids, "Raw gene IDs")
  counts <- SummarizedExperiment::assay(raw, "counts")
  ## prove the raw object reproduces the existing pseudobulks before altering spots.
  reconstruction <- list()
  for (context in names(inputs$objects)) {
    obj <- inputs$objects[[context]]
    expected <- obj$samples[obj$samples$SpD %in% cfg$spds, ]
    sel <- which(gm & labels[[context]] & keys %in% expected$key)
    groups <- match(keys[sel], expected$key)
    ix <- match(obj$genes$gene_id, gene_ids)
    assert(!anyNA(ix), "Original pseudobulk genes missing from raw object")
    membership <- Matrix::sparseMatrix(i = seq_along(sel), j = groups, x = 1,
                                       dims = c(length(sel), nrow(expected)))
    rebuilt_counts <- as.matrix(counts[ix, sel, drop = FALSE] %*% membership)
    expected_counts <- obj$counts[, expected$key, drop = FALSE]
    delta <- max(abs(rebuilt_counts - expected_counts))
    cells_agree <- identical(as.integer(tabulate(groups, nbins = nrow(expected))), as.integer(expected$ncells))
    reconstruction[[context]] <- data.frame(context = context, max_count_difference = delta,
      spot_counts_agree = cells_agree, samples = nrow(expected), genes = length(ix))
  }
  reconstruction <- do.call(rbind, reconstruction)
  write_tsv(reconstruction, file.path(env$out, "overlap", "raw_pseudobulk_reconciliation.tsv"))
  assert(all(reconstruction$max_count_difference < 1e-8 & reconstruction$spot_counts_agree),
         "Raw spot object does not reproduce existing pseudobulks; review input versions before reaggregation")
  statuses <- list()
  for (i in seq_len(nrow(screens))) {
    s <- screens[i, ]; hit_ids <- primary$gene_id[primary$screen_id == s$screen_id & primary$screen_hit]
    ## with overlap_refit_all, also track pairs whose only failed gate is the mediator gate.
    near_ids <- primary$gene_id[primary$screen_id == s$screen_id & primary$failed_gates %in% c("", "mediator_gate")]
    track_ids <- unique(c(hit_ids, if (isTRUE(cfg$overlap_refit_all)) near_ids))
    if (!length(hit_ids) && !isTRUE(cfg$overlap_refit_all)) {
      statuses[[i]] <- data.frame(screen_id = s$screen_id, status = "not_required_no_primary_hits", message = "", n_samples = NA_integer_, n_donors = NA_integer_)
      next
    }
    a <- get_artifact(s$screen_id, env)
    common <- labels[[s$source]] & labels[[s$target]]
    rebuilt <- inputs
    statuses[[i]] <- tryCatch({
      for (context in c(s$source, s$target)) {
        original <- inputs$objects[[context]]
        ix <- match(original$genes$gene_id, gene_ids)
        assert(!anyNA(ix), "Pseudobulk genes missing from raw object")
        sel <- which(gm & labels[[context]] & !common & keys %in% a$samples$key)
        groups <- match(keys[sel], a$samples$key)
        ncells <- tabulate(groups, nbins = nrow(a$samples))
        z <- Matrix::sparseMatrix(i = seq_along(sel), j = groups, x = 1,
                                  dims = c(length(sel), nrow(a$samples)))
        pb <- as.matrix(counts[ix, sel, drop = FALSE] %*% z)
        rownames(pb) <- original$genes$gene_id; colnames(pb) <- a$samples$key
        d <- a$samples; d$ncells <- ncells
        valid <- ncells >= cfg$min_spots & colSums(pb) > 0
        write_tsv(data.frame(key = d$key, ncells = ncells, retained = valid),
                  file.path(env$out, "overlap", paste0(s$screen_id, "_", context, "_retention.tsv")))
        pb <- pb[, valid, drop = FALSE]; d <- d[valid, ]
        assert(ncol(pb) > 0, "No pseudobulks remain after shared-spot removal")
        ## same normalization as the stored pseudobulk logcounts (spatialLIBD::registration_pseudobulk).
        lc <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(pb)), log = TRUE, prior.count = 1)
        rebuilt$objects[[context]] <- list(counts = pb, logcounts = lc, genes = original$genes, samples = d)
      }
      b <- run_screen(s, rebuilt, cfg, env, restrict = list(source = rownames(a$sdge), target = rownames(a$tdge)))
      r <- b$results; r$analysis <- "shared_spots_removed"; r$primary_hit <- r$gene_id %in% hit_ids
      write_tsv(r, file.path(env$out, "overlap", paste0(s$screen_id, "_results.tsv.gz")))
      cols <- c("gene_id", "screen_hit", "failed_gates", "a_beta", "a_p", "c_beta", "c_p", "cprime_beta", "cprime_p", "b_beta", "b_q", "relative_shrinkage")
      comparison <- merge(primary[primary$screen_id == s$screen_id & primary$gene_id %in% track_ids, c("gene_name", cols)],
        r[, cols], by = "gene_id", all.x = TRUE, suffixes = c("_primary", "_overlap_removed"))
      comparison$testable <- !is.na(comparison$screen_hit_overlap_removed)
      write_tsv(comparison, file.path(env$out, "overlap", paste0(s$screen_id, "_hit_comparison.tsv")))
      save_atomic(list(samples = b$samples, scaling = b$scaling, signature = env$signature), file.path(env$out, "overlap", paste0(s$screen_id, "_metadata.rds")))
      data.frame(screen_id = s$screen_id, status = "complete", message = "", n_samples = nrow(b$samples), n_donors = length(unique(b$samples$donor)))
    }, error = function(e) data.frame(screen_id = s$screen_id, status = "untestable_or_failed", message = conditionMessage(e), n_samples = NA_integer_, n_donors = NA_integer_))
  }
  status <- do.call(rbind, statuses)
  status$run_signature <- env$signature
  write_tsv(status, file.path(env$out, "overlap", "status.tsv"))
  raw_info <- data.frame(path = normalizePath(path), size = file.info(path)$size,
                         mtime = as.character(file.info(path)$mtime), md5 = unname(tools::md5sum(path)))
  write_tsv(raw_info, file.path(env$out, "overlap", "raw_input_fingerprint.tsv"))
}
