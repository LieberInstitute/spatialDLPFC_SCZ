#!/bin/env Rscript
library(here)
library(arrow)
library(data.table)
library(qs2)
library(stringi)

here::i_am('.git/HEAD')
refdir=here("processed-data", "ref")

fgwas <- file.path(refdir, 'GWAS_SCZ_tqtl-matched.tab.gz') ## rebuild if not found
if (file.exists(fgwas)) {
  gwas_tqtl <- fread(fgwas)
} else {
  stop("GWAS-tensorQTL variant-harmonized file not found: ", fgwas, "\nRun 03_eqtl_explore.Rmd to create it.")
}

## symlink the eqtl results from
eqtl_dir <- 'tqtl_out'
out_root <- here("processed-data", "eQTL", "coloc")

## prepare coloc input

library(coloc)
library(BiocParallel)

## inputs
## - gwas_tqtl: data.table with columns variant_id, BETA, SE, and NCAS/NCON (or s)
##    this must have GWAS SCZ variants harmonized to tensorQTL variant IDs
## - eqtl_dir: directory containing parquet files like:
##             PolyA.gene.chr1.parquet, PolyA.tx.chr1.parquet, RiboZ.gene.chr1.parquet, RiboZ.tx.chr1.parquet
## - out_root: where to save coloc RDS by window (w500k, w1mb)
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

## case fraction s for GWAS
if (all(c("ncas","ncon") %chin% names(gwas_tqtl))) {
  s_vec <- gwas_tqtl[, ncas/(ncas+ncon)]
  s <- stats::median(s_vec, na.rm = TRUE)
  stopifnot(stats::sd(s_vec, na.rm = TRUE) < 1e-2)
} else if ("s" %chin% names(gwas_tqtl)) {
  s <- unique(gwas_tqtl[!is.na(s), s])[1L]
} else {
  stop("Need ncas/ncon or a scalar column 's' in gwas_tqtl as_wide to set case fraction.")
}

## GWAS slice used for joins
stopifnot(all(c("variant_id","beta","beta_se") %chin% names(gwas_tqtl)))
gwas_idx <- gwas_tqtl[, .(snp = variant_id, beta_gwas = beta, varbeta_gwas = beta_se^2)]
setkey(gwas_idx, snp)

## helper: read all parquet for a prefix (e.g., "spd01.gene", "vasc.gene"), keep only needed cols
read_eqtl_nominal <- function(eqtl_dir, prefix, cis_window=500000) {
  files <- list.files(eqtl_dir, pattern = paste0("^", gsub("\\.", "\\\\.", prefix), ".*\\.parquet$"),
                      full.names = TRUE)
  if (length(files) == 0L) stop("No parquet files found for prefix: ", prefix)
  dt_list <- lapply(
    files,
    function(f) as.data.table(arrow::read_parquet(
      f,
      col_select = c("phenotype_id", "variant_id", "start_distance", "slope", "slope_se")
    ))
  )
  eqtl_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  ## filter to cis-window only
  eqtl_dt <- eqtl_dt[abs(start_distance) <= cis_window]
  ## ensure character to avoid factor surprises
  eqtl_dt[, `:=`(phenotype_id = as.character(phenotype_id),
                 variant_id   = as.character(variant_id))]
  eqtl_dt
}

## build coloc_df
build_coloc_df <- function(eqtl_dt, gwas_idx) {
  ## restrict to SNPs present in GWAS
  eqtl_dt <- eqtl_dt[variant_id %chin% gwas_idx$snp]
  ## rename/compute required fields
  eqtl_dt[, `:=`(
    snp = variant_id,
    beta_eqtl = slope,
    varbeta_eqtl = slope_se^2
  )]
  ## keep only the necessary columns
  eqtl_dt <- eqtl_dt[, .(phenotype_id, snp, beta_eqtl, varbeta_eqtl)]
  ## left join GWAS betas/vars
  setkey(eqtl_dt, snp)
  coloc_dt <- gwas_idx[eqtl_dt, on = "snp"]  ## joins beta_gwas/varbeta_gwas onto eqtl rows
  ## drop SNPs missing GWAS stats just in case
  coloc_dt <- coloc_dt[!is.na(beta_gwas) & !is.na(varbeta_gwas)]
  ## reorder columns
  setcolorder(coloc_dt, c("phenotype_id","snp","beta_eqtl","varbeta_eqtl","beta_gwas","varbeta_gwas"))
  coloc_dt[]
}

## harden run_coloc1: drop bad rows before ABF, keep priors visible
run_coloc1 <- function(gene, coloc_dt, s, sdY = 1,
                       p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, min_snps = 10L) {
  d <- coloc_dt[phenotype_id == gene]
  ## drop NA/Inf and zero-variance rows up front
  d <- d[is.finite(beta_eqtl) & is.finite(varbeta_eqtl) & varbeta_eqtl > 0 &
         is.finite(beta_gwas) & is.finite(varbeta_gwas) & varbeta_gwas > 0]
  if (nrow(d) < min_snps) return(NULL)
  ## optional quick signal check
  z_ok <- any(abs(d$beta_eqtl / sqrt(d$varbeta_eqtl)) >= 2, na.rm = TRUE)
  if (!z_ok) return(NULL)
  if (any(duplicated(d$snp))) d <- unique(d, by = "snp")
  eqtl <- list(snp = d$snp, beta = d$beta_eqtl, varbeta = d$varbeta_eqtl,
               type = "quant", sdY = sdY)
  gwas <- list(snp = d$snp, beta = d$beta_gwas, varbeta = d$varbeta_gwas,
               type = "cc", s = s)
  coloc::check_dataset(eqtl)
  coloc::check_dataset(gwas)
  coloc::coloc.abf(eqtl, gwas, p1 = p1, p2 = p2, p12 = p12)
}

## wrap with tryCatch so one bad gene doesn't kill the batch
run_coloc1_safe <- function(g, coloc_dt, s, ...) {
#  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table not available on worker")
#  if (!requireNamespace("coloc", quietly = TRUE)) stop("coloc not available on worker")
#  library(data.table)
#  library(coloc)
  tryCatch(run_coloc1(g, coloc_dt = coloc_dt, s = s, ...),
           error = function(e) structure(list(.error = conditionMessage(e)), class = "coloc_err"))
}

n_cores = 8
bpparams <- MulticoreParam(n_cores, stop.on.error = FALSE, progressbar = FALSE)
#BiocParallel::bpstart(bpparams)
#on.exit(BiocParallel::bpstop(param_ps), add = TRUE)

runs <- paste0(c(paste0('spd0',1:7), "neun", "vasc", "pnn", "neuropil"), ".gene")
window_cfg <- c(w500k = 500000L, w1mb = 1000000L)

priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, sdY = 1, min_snps = 10L)

for (wtag in names(window_cfg)) {
  cis_window <- as.integer(window_cfg[[wtag]])
  out_dir <- file.path(out_root, wtag)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  message(Sys.time(), " | window=", wtag, " | cis_window=", cis_window)

  for (prefix in runs) {
    t0 <- Sys.time()
    ## determine output path first so we can skip completed runs
    out_path <- file.path(out_dir, paste0("coloc_", gsub("\\.", "_", prefix), ".qs2"))
    meta_path <- sub("\\.qs2$", ".runmeta.tsv.gz", out_path)
    if (file.exists(out_path) && file.exists(meta_path)) {
      message(Sys.time(), " | ", wtag, " | ", prefix, " | output+meta exist, skipping")
      next
    }

    message(Sys.time(), " | ", wtag, " | ", prefix, " | reading eQTL parquet...")
    eqtl_dt <- read_eqtl_nominal(eqtl_dir, prefix, cis_window = cis_window)
    n_eqtl_in_window <- nrow(eqtl_dt)

    message(Sys.time(), " | ", wtag, " | ", prefix, " | building coloc_df...")
    coloc_dt <- build_coloc_df(eqtl_dt, gwas_idx)
    n_rows_coloc <- nrow(coloc_dt)
    n_snps_coloc <- uniqueN(coloc_dt$snp)

    genes <- unique(as.character(coloc_dt$phenotype_id))
    genes <- genes[!is.na(genes)]
    n_genes_input <- length(genes)
    message(Sys.time(), " | ", wtag, " | ", prefix, " | running coloc per feature: ", n_genes_input, " genes")

    res_list <- bplapply(setNames(genes, genes), run_coloc1_safe,
      coloc_dt = coloc_dt,
      s = s,
      p1 = priors$p1,
      p2 = priors$p2,
      p12 = priors$p12,
      sdY = priors$sdY,
      min_snps = priors$min_snps,
      BPPARAM = bpparams
    )

    ## drop errors and NULLs before naming/saving
    err_idx <- which(vapply(res_list, inherits, logical(1), "coloc_err"))
    n_err <- length(err_idx)
    if (n_err) {
      msg <- unique(vapply(res_list[err_idx], `[[`, character(1), ".error"))
      message(n_err, " genes threw coloc errors; first: ", msg[1])
      res_list <- res_list[-err_idx]
    }
    null_idx <- vapply(res_list, is.null, logical(1))
    n_null <- sum(null_idx)
    res_list <- res_list[!null_idx]
    n_saved <- length(res_list)

    qs_save(res_list, out_path)
    elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    meta <- data.table(
      timestamp = as.character(Sys.time()),
      context = prefix,
      window_tag = wtag,
      cis_window = cis_window,
      out_path = out_path,
      n_eqtl_in_window = n_eqtl_in_window,
      n_rows_coloc = n_rows_coloc,
      n_snps_coloc = n_snps_coloc,
      n_genes_input = n_genes_input,
      n_genes_error = n_err,
      n_genes_null = n_null,
      n_genes_saved = n_saved,
      priors_p1 = priors$p1,
      priors_p2 = priors$p2,
      priors_p12 = priors$p12,
      priors_sdY = priors$sdY,
      min_snps = priors$min_snps,
      case_fraction_s = s,
      elapsed_sec = elapsed_sec
    )
    fwrite(meta, meta_path, sep = "\t")
    message(Sys.time(), " | ", wtag, " | ", prefix, " | saved: ", out_path)
  }
}
