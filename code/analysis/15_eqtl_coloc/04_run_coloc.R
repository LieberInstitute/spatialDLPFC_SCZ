#!/bin/env Rscript
library(here)
library(arrow)
library(data.table)
library(bigsnpr) # remotes::install_github("privefl/bigsnpr")
library(bigreadr) ## for fread2
library(qs2)
library(stringi)

here::i_am('.git/HEAD')
refdir=here("processed-data", "ref")

## tensorqtl variant universe:

pvar <- fread(here("processed-data/genotypes/plink2", "merged_maf05.pvar"), sep = "\t")
setnames(pvar, '#CHROM', 'CHROM')

## symlink the eqtl results from 06_eQTL
eqtl_dir <- 'tqtl_out'
out_dir  <- 'coloc'
info_snp <- pvar[, .(chr = as.character(CHROM),
  pos = as.integer(POS),
  a0  = REF, a1  = ALT, ID)]

fgwas <- file.path(refdir, 'scz_gwas_wide.scz.tqtl-matched.tab.gz') ## rebuild if not found
liftOver <- Sys.which("liftOver")
stopifnot(file.exists(liftOver))
if (file.exists(fgwas)) {
  gwas_wide <- fread(fgwas)
} else {
  fhg19gwas <- file.path(refdir, 'PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz')
  #data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz'
  stopifnot(file.exists(fhg19gwas))
  gwas_wide <- fread2(fhg19gwas)
  colnames(gwas_wide) <- tolower(colnames(gwas_wide))
  setnames(gwas_wide, old = c("chrom","id","a1","a2","se",    "neff","pval"),
                      new = c("chr","rsid","a0","a1","beta_se","N",  "p"))
  setcolorder(gwas_wide, c("rsid","chr","pos","a0","a1","beta","beta_se","N","p"))
  gwas_wide$chr <- as.character(gwas_wide$chr)

  gwas_wide <- snp_modifyBuild(gwas_wide, liftOver, from = 'hg19', to = 'hg38')
  ## 4967 variants have not been mapped.
  library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
  dbsnps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
  setDT(gwas_wide)
  gwas_wide$chr <- as.character(gwas_wide$chr)
  ## only pull the unmapped unique rsids:
  ids_unmapped <- unique(gwas_wide[is.na(pos) & grepl("^rs\\d+$", rsid), rsid])
  ## annoyingly slow:
  rehg38 <- snpsById(dbsnps, ids=ids_unmapped, ifnotfound="drop")
  ## rescued 4360 variants
  ## Normalize the SNPlocs result to chr/pos and join-update only rows with missing pos
  hdt <- as.data.table(rehg38)
  setnames(hdt, c('seqnames'), c('chr'))
  hdt$chr <- as.character(hdt$chr)
  setkey(hdt, RefSNP_id)
  setkey(gwas_wide, rsid)
  gwas_wide[hdt$RefSNP_id, `:=`(chr=hdt$chr, pos=hdt$pos)]
   ## only 607 variants are not mapped
  gwas_wide[, variant_id := sprintf("chr%s:%s:%s:%s", chr, pos, a0, a1)]
  gwas_wide[, chr := paste0("chr", chr)]
  ## harmonize with tensorqtl variant universe
  ## prepare GWAS for matching
  sumstats_for_match <- gwas_wide[, .(
    chr, pos,
    a0 = a1,      # A2
    a1 = a0,      # A1 (effect allele)
    beta, beta_se, N, p, ncas, ncon, impinfo, fcas, fcon, rsid
  )]

  gwas_m <- bigsnpr::snp_match(sumstats_for_match, info_snp, return_flip_and_rev = TRUE)
  # gwas_m now has alleles aligned to info_snp, and beta flipped if alleles swapped
  setDT(gwas_m)
  # Rebuild variant_id to match tensorQTL
  #gwas_m[, variant_id := sprintf("%s:%d:%s:%s", chr, pos, a0, a1)]
  gwas_m[, variant_id := info_snp$ID[`_NUM_ID_`]]
  ## build the join key
  gwas_fixed <- gwas_m[, .(
    variant_id, beta, beta_se,
    N,  p,  ncas, ncon, impinfo,
    rsid, fcas, fcon
  )]
  # deduplicate defensively
  setkey(gwas_fixed, variant_id)
  gwas_fixed <- unique(gwas_fixed)
  ##
  fwrite(gwas_fixed, fgwas, sep = "\t")
}

## prepare coloc input

library(coloc)
library(BiocParallel)

## inputs
## - gwas_wide: data.table with columns variant_id, BETA, SE, and NCAS/NCON (or s)
## - eqtl_dir: directory containing parquet files like:
##             PolyA.gene.chr1.parquet, PolyA.tx.chr1.parquet, RiboZ.gene.chr1.parquet, RiboZ.tx.chr1.parquet
## - out_dir: where to save coloc RDS per run
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## case fraction s for GWAS
if (all(c("ncas","ncon") %chin% names(gwas_wide))) {
  s_vec <- gwas_wide[, ncas/(ncas+ncon)]
  s <- stats::median(s_vec, na.rm = TRUE)
  stopifnot(stats::sd(s_vec, na.rm = TRUE) < 1e-2)
} else if ("s" %chin% names(gwas_wide)) {
  s <- unique(gwas_wide[!is.na(s), s])[1L]
} else {
  stop("Need ncas/ncon or a scalar column 's' in gwas_wide to set case fraction.")
}

## GWAS slice used for joins
stopifnot(all(c("variant_id","beta","beta_se") %chin% names(gwas_wide)))
gwas_idx <- gwas_wide[, .(snp = variant_id, beta_gwas = beta, varbeta_gwas = beta_se^2)]
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

n_cores = 6
bpparams <- MulticoreParam(n_cores, stop.on.error = FALSE, progressbar = FALSE)
#BiocParallel::bpstart(bpparams)
#on.exit(BiocParallel::bpstop(param_ps), add = TRUE)

runs <- paste0(c(paste0('spd0',1:7), "neun", "vasc", "pnn", "neuropil"), ".gene")

for (prefix in runs) {
  ## determine output path first so we can skip completed runs
  out_path <- file.path(out_dir, paste0("coloc_", gsub("\\.", "_", prefix), ".qs2"))
  if (file.exists(out_path)) {
    message(Sys.time(), " | ", prefix, " | output exists, skipping: ", out_path)
    next
  }
  message(Sys.time(), " | ", prefix, " | reading eQTL parquet...")
  eqtl_dt <- read_eqtl_nominal(eqtl_dir, prefix)
  message(Sys.time(), " | ", prefix, " | building coloc_df...")
  coloc_dt <- build_coloc_df(eqtl_dt, gwas_idx)
  message(Sys.time(), " | ", prefix, " | running coloc per feature...")

  genes <- unique(as.character(coloc_dt$phenotype_id))
  genes <- genes[!is.na(genes)]

  res_list <- bplapply( setNames(genes, genes), run_coloc1_safe,
    coloc_dt = coloc_dt,  s = s,  BPPARAM = bpparams
  )

  ## drop errors and NULLs before naming/saving
  err_idx <- which(vapply(res_list, inherits, logical(1), "coloc_err"))
  if (length(err_idx)) {
    msg <- unique(vapply(res_list[err_idx], `[[`, character(1), ".error"))
    message(length(err_idx), " genes threw coloc errors; first: ", msg[1])
    res_list <- res_list[-err_idx]
  }
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]

  qs_save(res_list, out_path)
  message(Sys.time(), " | ", prefix, " | saved: ", out_path)
}
