#!/bin/env Rscript

library("SpatialExperiment")
library("here")
here::i_am('.git/HEAD')
datadir=here("processed-data", "rds")
library(sessioninfo)
library(tidyverse)
library(data.table)
library(qs2)
library(sva)
library(edgeR)
library(matrixStats)

## strategy: dsarg must be one of "spd" or "spg"
## use optargs to take this as the first argument of the script
dsarg <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.null(dsarg) && nchar(dsarg)>1)

modelstr='~DX + sex + age + slide_id + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5'

if (tolower(dsarg) %in% c('spd', 'spg')) {
  dsarg <- tolower(dsarg)
} else {
  stop("First argument must be one of 'spd' or 'spg'")
}

if (dsarg=='spd') {
  message("Preparing eQTL input files for spatial domains")
  d <- readRDS(file.path(datadir, "07_dx_pseudobulk", "sce_pseudo_PRECAST07_donor_spd.rds"))
  spd_indices <- split(seq_len(ncol(d)), colData(d)$registration_variable)
  dlist <- lapply(spd_indices, function(idx) d[, idx])
} else  { # mist be spg
  message("Preparing eQTL input files for spatial proteogenomics microenvironments")
  flist <- list(neuropil="pseudo_neuropil_pos_donor_spd.rds",
                neun="pseudo_neun_pos_donor_spd.rds",
                pnn="pseudo_pnn_pos_donor_spd.rds"
                # vasc="pseudo_vasc_pos_donor_spd.rds" -- missing slide_id
              )
  dlist <- lapply(flist, function(x) {
      spe <- readRDS(file.path(datadir, "PB_dx_spg", x))
      agg <- scuttle::aggregateAcrossCells(spe, DataFrame(donor = spe$brnum), store_number = "n_SpDs")
      ## colnames are lost during aggregation, restore them from donor
      ncells_sum <- rowsum(
        x = matrix(spe$ncells, ncol = 1),
        group = spe$brnum,  na.rm = TRUE
      )
      donor_ids <- colData(agg)$donor
      colData(agg)$ncells <- ncells_sum[donor_ids, 1]
      agg
  })
}

dlist <- lapply(dlist, function(spe) {
  ## make sure brnum values are unique
  stopifnot(length(unique(spe$brnum))==ncol(spe))
  colnames(spe) <- spe$brnum
  return(spe)
})


snpPCs <- NULL
## path to SNP PCs (eigenvec) prepared with get_SNP_PCs.sh
fsnp_pcs <- 'genotypes/plink2/merged_maf05_pca.eigenvec'
stopifnot(file.exists(fsnp_pcs)) ## the RSEs must be symlinked or copied in data/
if (!is.null(fsnp_pcs)) {
  snpPCs <- fread(fsnp_pcs, data.table=F)
  if (is.null(snpPCs$SAMPLE_ID)) {
    ## must determine if FID is also present, not just #IID
    nidcols <- ncol(snpPCs) - length(grep('PC\\d+$', colnames(snpPCs)))
    if (nidcols==2) snpPCs <- snpPCs[ , -1] ## FID is first, discard the first column
    colnames(snpPCs)[1]  <- 'SAMPLE_ID'
  }
  rownames(snpPCs) <- snpPCs$SAMPLE_ID ## this is the genotype ID
  colnames(snpPCs) <- gsub('^PC', 'snpPC', colnames(snpPCs))
}


prep_spe_qtl <- function(
  spe,  outlier_mad = 0,  use_ncells = TRUE, ## disable qc filtering by setting outlier_mad=0
  ## BED gene detectability filter (scRNA-style 20%)
  bed_detect_prop = 0.20,
  bed_min_cpm = 0.10,
  bed_min_count = 6,
  min_detect_n_floor = 3,
  ## PCA gene filter (paper-style pi0)
  pca_pi0_max = 0.90,
  pca_use_hvg = TRUE,
  pca_use_detect_filter = TRUE,
  pca_hvg_n = 5000, ## use top 5000 HVGs for PCA
  ## transformations
  prior_count = 1,
  ## minimum sample size checks
  min_n_qtl_hard = 20,
  min_resid_df = 20,
  covar_formula_known = NULL,  # e.g. ~ Dx + Age + Sex + slideID + snpPC1 + ... + snpPC5
  stop_if_too_small = FALSE,
  verbose = TRUE
) {
  stopifnot("counts" %in% assayNames(spe))
  cnt0 <- assay(spe, "counts")

  rep <- data.frame(
    step = character(),
    n_genes = integer(),
    n_samples = integer(),
    removed_genes = integer(),
    removed_samples = integer(),
    stringsAsFactors = FALSE
  )

  add_rep <- function(step, spe_obj, rg = 0L, rs = 0L) {
    rep <<- rbind(rep, data.frame(
      step = step,
      n_genes = nrow(spe_obj),
      n_samples = ncol(spe_obj),
      removed_genes = as.integer(rg),
      removed_samples = as.integer(rs),
      stringsAsFactors = FALSE
    ))
    if (verbose) {
      message(sprintf(
        "%-28s genes=%6d  samples=%4d  (-genes=%4d, -samples=%3d)",
        step, nrow(spe_obj), ncol(spe_obj), rg, rs
      ))
    }
  }
  add_rep("start", spe)
  ## Sample QC filtering?
  if (outlier_mad > 0) {
    lib <- colSums(cnt0)
    lib_z <- (lib - median(lib)) / mad(lib, constant = 1)
    keep <- lib_z > (-outlier_mad)

    if (use_ncells && ("ncells" %in% colnames(colData(spe)))) {
      nc <- spe$ncells
      if (!all(is.na(nc))) {
        nc_z <- (nc - median(nc, na.rm = TRUE)) / mad(nc, constant = 1, na.rm = TRUE)
        keep <- keep & (nc_z > (-outlier_mad))
      }
    }
    spe_qc <- spe[, keep, drop = FALSE]
    add_rep("sample_qc_mad", spe_qc, rg = 0L, rs = sum(!keep))
  } else {
    spe_qc <- spe
  }
  n <- ncol(spe_qc)
  if (n < min_n_qtl_hard) {
    msg <- sprintf("n=%d < hard minimum %d donors for QTL mapping (very low power).", n, min_n_qtl_hard)
    if (stop_if_too_small) stop(msg) else if (verbose) message("WARNING: ", msg)
  }
  ## optional: DF-based feasibility check for known covariates
  known_rank <- NA_integer_
  resid_df_known <- NA_integer_
  if (!is.null(covar_formula_known)) {
    mm <- model.matrix(covar_formula_known, data = as.data.frame(colData(spe_qc)))
    known_rank <- qr(mm)$rank
    resid_df_known <- n - known_rank - 1L
    if (verbose) {
      message(sprintf("known covariates: rank=%d, residual_df (before exprPCs)=%d", known_rank, resid_df_known))
    }
    if (!is.na(resid_df_known) && resid_df_known < min_resid_df) {
      msg <- sprintf("Residual DF too small with known covariates: %d < %d. Reduce covariates / merge levels / increase n.",
                     resid_df_known, min_resid_df)
      if (stop_if_too_small) stop(msg) else if (verbose) message("WARNING: ", msg)
    }
  }
  ## BED gene filter (20% detectable)
  cnt <- assay(spe_qc, "counts")
  dge0 <- DGEList(cnt)
  dge0 <- calcNormFactors(dge0, method = "TMM")
  cpm0 <- cpm(dge0, log = FALSE)

  min_detect_n <- max(ceiling(bed_detect_prop * n), min_detect_n_floor)

  ## detectable donors per gene under both CPM and raw count thresholds
  det <- (cpm0 >= bed_min_cpm) & (cnt >= bed_min_count)
  det_n <- rowSums(det)

  keep_gene_bed <- det_n >= min_detect_n
  spe_bed <- spe_qc[keep_gene_bed, , drop = FALSE]
  add_rep(sprintf("bed_gene_filter_%dpct", as.integer(100 * bed_detect_prop)),
          spe_bed, rg = sum(!keep_gene_bed), rs = 0L)
  ## BED normalization + INT (tensorQTL phenotype)
  dge_bed <- DGEList(assay(spe_bed, "counts"))
  dge_bed <- calcNormFactors(dge_bed, method = "TMM")
  log2cpm <- cpm(dge_bed, log = TRUE, prior.count = prior_count)

  rr <- rowRanks(log2cpm, ties.method = "average")
  pp <- (rr - 0.5) / ncol(log2cpm)
  log2cpm_rint <- qnorm(pp)

  assays(spe_bed)$tmm <- log2cpm
  assays(spe_bed)$rint <- log2cpm_rint
  add_rep("bed_TMM_log2cpm_INT", spe_bed, rg = 0L, rs = 0L)

  ## PCA gene filter (pi0)
  cnt_qc <- assay(spe_qc, "counts")
  pi0 <- rowMeans(cnt_qc == 0)
  if (verbose) {
    qs <- quantile(pi0, probs = c(0, 0.5, 0.9, 0.95, 0.99, 1))
    message(sprintf("pi0 quantiles: min=%.3f med=%.3f p90=%.3f p95=%.3f p99=%.3f max=%.3f",
                    qs[1], qs[2], qs[3], qs[4], qs[5], qs[6]))
  }
  keep_gene_pca <- pi0 < pca_pi0_max
  if (pca_use_detect_filter) {
    keep_gene_pca <- keep_gene_pca & (det_n >= min_detect_n)
  }
  spe_pca <- spe_qc[keep_gene_pca, , drop = FALSE]
  add_rep(sprintf("pca_gene_filter_pi0<%.2f%s",
                pca_pi0_max, if (pca_use_detect_filter) "_and_detect" else ""),
          spe_pca, rg = sum(!keep_gene_pca), rs = 0L)
  ## PCA matrix: TMM + log2CPM, then z-score per gene
  dge_pca <- DGEList(assay(spe_pca, "counts"))
  dge_pca <- calcNormFactors(dge_pca, method = "TMM")
  log2cpm_pca <- cpm(dge_pca, log = TRUE, prior.count = prior_count)
  if (pca_use_hvg) {
    cpm_lin <- cpm(dge_pca, log = FALSE)
    mu <- rowMeans(cpm_lin)
    va <- rowVars(cpm_lin)
    fano <- va / pmax(mu, 1e-8)
    o <- order(fano, decreasing = TRUE)
    sel <- o[seq_len(min(pca_hvg_n, length(o)))]
    log2cpm_pca <- log2cpm_pca[sel, , drop = FALSE]
    prev_g <- nrow(spe_pca)
    spe_pca <- spe_pca[sel, , drop = FALSE]
    add_rep(sprintf("pca_HVG_top%d", pca_hvg_n), spe_pca, rg = prev_g - nrow(spe_pca), rs = 0L)
  }

  mat_z <- t(scale(t(log2cpm_pca), center = TRUE, scale = TRUE))
  mat_z[!is.finite(mat_z)] <- 0
  assays(spe_pca)$tmm <- log2cpm_pca
  assays(spe_pca)$zscore <- mat_z
  #pc <- prcomp(t(mat_z), center = FALSE, scale. = FALSE)
  #expr_pcs <- pc$x
  #colnames(expr_pcs) <- paste0("exprPC", seq_len(ncol(expr_pcs)))

  ## overall feasibility summary
  qtl_ok <- (n >= min_n_qtl_hard)
  list(
    spe_bed = spe_bed,
    spe_pca = spe_pca,
    #expr_pcs = expr_pcs,
    report = rep,
    settings = list(
      n = n,
      min_detect_n = min_detect_n,
      bed_detect_prop = bed_detect_prop,
      bed_min_cpm = bed_min_cpm,
      bed_min_count = bed_min_count,
      pca_pi0_max = pca_pi0_max,
      known_rank = known_rank,
      resid_df_known = resid_df_known
    ),
    qtl_ok = qtl_ok
  )
}


## make the
spe2bed <- function(spe_bed) { ## assumes tmm assay is present
  rr_df <- as.data.frame(rowRanges(spe_bed))
  tmm <- assays(spe_bed)$tmm ## already TMM normalized + INT
  stopifnot(!is.null(tmm))
  ## best practice for tensorQTL is to use rank-inverse normal transformation (RINT)
  rinvnorm <- function(x) qnorm((rank(x, ties.method='average') - 0.5) / length(x))
  counts <- t(apply(tmm, 1, rinvnorm)) ## BED uses samples in rows
  #counts <- t(rinvnorm)
  ## prepare BED data frame
  rr_df$tss_start=ifelse(rr_df$strand=='-', rr_df$end, rr_df$start)
  rr_df$start <- rr_df$tss_start
  rr_df <- rr_df %>%  tibble::rownames_to_column("ID") %>%
  dplyr::arrange(seqnames, start) %>%
  dplyr::mutate(end = start + 1) %>%
  dplyr::select(`#Chr` = seqnames, start, end, ID)
  colnames(counts) <- colnames(spe_bed) ## brnum
  counts <- as.data.frame(counts) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(rr_df, ., by = "ID")
  return(counts)
}

covar_format <- function(data) {
  data <- as.data.frame(data)
  data <- t(data)
  data <- as.data.frame(data) %>% rownames_to_column("id")
  return(data)
}


odir='tqtl_in' ## output directory for expression and covariate files

if (!dir.exists(odir)) dir.create(odir)

model <- NULL

## optional: reduce frses to only a subset of features
#frses <- frses[c('gene')]

for (i in seq_len(length(dlist))) {

  ## start filtering and transforming the data
  dsname <-names(dlist)[[i]]
  spe <- dlist[[i]]

  fn  <- paste0(odir, '/', dsname, '.gene.exprPCs.qs2')
  fncovars <- sub('exprPCs.qs2', 'covars.txt', fn, fixed = TRUE)
  fbed  <- sub('exprPCs.qs2', 'expr.bed.gz', fn, fixed = TRUE)
  ## if all these files exist, skip
  if (file.exists(fn) && file.exists(fncovars) && file.exists(fbed)) {
    cat(".. all files for", dsname, "already exist, skipping\n")
    next
  }
  cat("Processing dataset", dsname, " ..\n")
  ## add snpPCs to colData, matching by brnum / SAMPLE_ID
  stopifnot(!is.null(snpPCs))
  stopifnot(all(spe$brnum %in% rownames(snpPCs)))
  colData(spe) <- cbind(colData(spe), snpPCs[spe$brnum, -1, drop = FALSE])

  pd <- as.data.frame(colData(spe))
  if (is.null(pd$DX)) {
    pd$DX <- toupper(pd$dx)
  }
  pd$DX <- factor(as.character(pd$DX), levels = c('NTC', 'SCZ'))

  ## model (sample) factors and covariates
  model <- model.matrix(as.formula(modelstr), data = pd)
  #mod0 <- model.matrix(as.formula(mod0str), data = pd)
  stopifnot(identical(rownames(model), rownames(pd)))
  stopifnot(identical(colnames(spe), rownames(pd)))
  ## discard chrY and chrM genes
  spe <- spe[!grepl('chr[YM]', seqnames(rowRanges(spe))), ]
  ## filter spe and produce spe_bed and spe_pca and expr_pcs

  res <- prep_spe_qtl(
    spe,
    covar_formula_known = as.formula(modelstr),
    verbose = TRUE
  )
  # res$report
  #expr_for_bed <- assays(res$spe_bed)$rint
  spe_bed <- res$spe_bed ## to use for BED export
  spe_pca <- res$spe_pca
  #expr_pcs <- res$expr_pcs

  ## calculate feature PCs (sva)
  ffPCs <- NULL
  ## note: set this to FALSE if we want to use log(rpkm+1) instead of existing logcounts
  #have_logcounts <- !is.null(assays(rse)$logcounts)
  have_zscore <- !is.null(assays(spe_pca)$zscore)
  stopifnot(have_zscore)
  if (file.exists(fn)) {
    cat("Loading", dsname, "gene expression PCs from", fn, "\n")
    ffPCs <- qs_read(fn)
  } else {
    cat("Calculate", dsname, "gene expression PCs..\n")
    lmx <- assays(spe_pca)$zscore
    pca <- prcomp(t(lmx))
    #k <- 0
    k <- num.sv(lmx, model)
    cat("  ..using", k, "expression PCs\n")
    ffPCs <- pca$x[, 1:k]
    #svobj <- sva::sva(lmx, mod=model, mod0=mod0, n.sv=k)
    #ffPCs<- svobj$sv
    colnames(ffPCs) <- paste0('exprPC', seq_len(ncol(ffPCs)))
    rownames(ffPCs) <- colnames(lmx)
    qs_save(ffPCs, file=fn)
    cat("     saved as",fn, "\n")
  }

  stopifnot(identical(rownames(ffPCs), colnames(spe_pca)))

  cat(" exporting data for",dsname,"..\n")

  cpd <- covar_format(model)
  ## PC data, expr PCs covariates
  stopifnot(identical(rownames(ffPCs), rownames(pd)))
  cpc <- covar_format(ffPCs)
  ## bind model (sample) covars with feature covars (SVs) and save
  stopifnot(identical(colnames(cpd), colnames(cpc)))
  covars <- rbind(cpd, cpc)
  #fn=paste0(odir, '/', dsname, '.', feature,'.featurePCs.qs')

  fwrite(covars, file = fncovars,
         sep = "\t", quote = FALSE, row.names = FALSE)
  bed <- spe2bed(spe_bed)
  fwrite(bed, file = fbed, sep = "\t", quote = FALSE, row.names = FALSE)
  rm(bed)
  cat("  ..exported to", fbed, "\n")
} ## for each spe object
