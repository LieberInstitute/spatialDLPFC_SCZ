suppressMessages(library(data.table))

default_eqtl_contexts <- c(sprintf("spd%02d", 1:7), "neun", "vasc", "pnn", "neuropil")

assert_file <- function(path) {
  if (!file.exists(path)) stop("missing file: ", path)
  invisible(path)
}

assert_files_exist <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("Missing required files:\n", paste(missing, collapse = "\n"))
  }
  invisible(files)
}

assert_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0L) {
    stop(sprintf("%s is missing required columns: %s", label, paste(miss, collapse = ", ")))
  }
  invisible(dt)
}

assert_cols_exist <- assert_cols

normalize_gene_id <- function(x) sub("\\.\\d+$", "", as.character(x))
norm_gene_id <- normalize_gene_id

drop_pandas_index_cols <- function(dt) {
  if ("V1" %in% names(dt)) dt[, V1 := NULL]
  unnamed <- grep("^Unnamed", names(dt), value = TRUE)
  if (length(unnamed) > 0L) dt[, (unnamed) := NULL]
  dt
}

warn_missing_rsid <- function(dt, label, rsid_col = "rsID") {
  nmiss <- dt[is.na(get(rsid_col)) | get(rsid_col) == "", .N]
  if (nmiss > 0L) {
    warning(sprintf("%s has %d rows with no %s annotation; keeping missing values.", label, nmiss, rsid_col))
  }
  invisible(dt)
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x)) x[[1]] else NA_character_
}

uniq_concat <- function(x) {
  x <- as.character(x)
  x <- sort(unique(x[!is.na(x) & x != ""]))
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

std_symbol_key <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y == ""] <- NA_character_
  y
}

get_deg_file_paths <- function(project_root = here::here(), rdsdir = file.path(project_root, "processed-data", "rds")) {
  list(
    spd = file.path(rdsdir, "10_dx_deg_adjust_spd", "dx-deg_PRECAST07.csv"),
    spg = list(
      neun = file.path(project_root, "code", "analysis", "dx_deg_spg_neun", "neun-dx_DEG-GM.csv"),
      neuropil = file.path(project_root, "code", "analysis", "dx_deg_spg_neuropil", "neuropil-dx_DEG-GM.csv"),
      pnn = file.path(project_root, "code", "analysis", "dx_deg_spg_pnn", "pnn-dx_DEG-GM.csv"),
      vasc = file.path(project_root, "processed-data", "spg_pb_de", "test_SPD_pseudo_vasc_pos.csv")
    )
  )
}

load_deg_sets <- function(
  contexts = default_eqtl_contexts,
  project_root = here::here(),
  rdsdir = file.path(project_root, "processed-data", "rds"),
  normalize_ids = TRUE
) {
  deg_paths <- get_deg_file_paths(project_root = project_root, rdsdir = rdsdir)
  assert_file(deg_paths$spd)
  assert_files_exist(unlist(deg_paths$spg, use.names = FALSE))

  id_fun <- if (normalize_ids) normalize_gene_id else as.character
  strict <- list()
  loose <- list()

  spd_degs <- fread(deg_paths$spd)
  assert_cols(spd_degs, c("ensembl", "fdr_scz", "p_value_scz"), basename(deg_paths$spd))
  spd_strict <- unique(id_fun(spd_degs[fdr_scz <= 0.1 & !is.na(ensembl), ensembl]))
  spd_loose <- unique(id_fun(spd_degs[p_value_scz < 0.05 & !is.na(ensembl), ensembl]))
  for (ctx in grep("^spd", contexts, value = TRUE)) {
    strict[[ctx]] <- spd_strict
    loose[[ctx]] <- spd_loose
  }

  for (ctx in names(deg_paths$spg)) {
    ddt <- fread(deg_paths$spg[[ctx]])
    assert_cols(ddt, c("ensembl", "fdr_scz", "p_value_scz"), basename(deg_paths$spg[[ctx]]))
    strict[[ctx]] <- unique(id_fun(ddt[fdr_scz <= 0.1 & !is.na(ensembl), ensembl]))
    loose[[ctx]] <- unique(id_fun(ddt[p_value_scz < 0.05 & !is.na(ensembl), ensembl]))
  }

  if (!setequal(names(strict), contexts) || !setequal(names(loose), contexts)) {
    stop("DGE context mismatch. Strict: ", paste(sort(names(strict)), collapse = ", "),
         " | Loose: ", paste(sort(names(loose)), collapse = ", "))
  }

  list(strict = strict, loose = loose, paths = deg_paths)
}

apply_deg_flags <- function(dt, context_col, gene_col, deg_sets, strict_col = "DGE", loose_col = "DGEp05") {
  d <- copy(dt)
  d[, .gene_norm_tmp__ := normalize_gene_id(get(gene_col))]
  for (ctx in unique(as.character(d[[context_col]]))) {
    d[get(context_col) == ctx, (c(strict_col, loose_col)) := list(
      as.integer(.gene_norm_tmp__ %chin% deg_sets$strict[[ctx]]),
      as.integer(.gene_norm_tmp__ %chin% deg_sets$loose[[ctx]])
    )]
  }
  d[, .gene_norm_tmp__ := NULL]
  d
}

apply_gwas_flag <- function(dt, variant_col, gwas_variants, out_col = "GWAS") {
  d <- copy(dt)
  d[, (out_col) := as.integer(as.character(get(variant_col)) %chin% as.character(gwas_variants))]
  d
}

build_strong_coloc_lookup <- function(coloc_res, context_col = "ds", gene_col = "fid", variant_col = "lead_snp") {
  assert_cols(coloc_res, c(context_col, gene_col, variant_col), "coloc results")
  if ("cat" %in% names(coloc_res)) {
    strong_coloc <- coloc_res[cat == "strong_coloc"]
  } else {
    assert_cols(coloc_res, c("PP4"), "coloc results")
    warning("coloc results have no cat column; using PP4 > 0.8 fallback.")
    strong_coloc <- coloc_res[PP4 > 0.8]
  }
  list(
    gene = unique(strong_coloc[, .(
      context = as.character(get(context_col)),
      gene_id_norm = normalize_gene_id(get(gene_col))
    )]),
    variant = unique(strong_coloc[, .(
      context = as.character(get(context_col)),
      variant_id = as.character(get(variant_col))
    )])
  )
}

apply_strong_coloc_flags <- function(
  dt,
  context_col,
  gene_col,
  variant_col,
  strong_lookup,
  gene_flag_col = "gene_coloc",
  variant_flag_col = "var_coloc"
) {
  d <- copy(dt)
  d[, .gene_norm_tmp__ := normalize_gene_id(get(gene_col))]
  d[, (variant_flag_col) := 0L]
  d[, (gene_flag_col) := 0L]
  if (nrow(strong_lookup$variant) > 0L) {
    d[strong_lookup$variant, on = setNames(c("context", "variant_id"), c(context_col, variant_col)), (variant_flag_col) := 1L]
  }
  if (nrow(strong_lookup$gene) > 0L) {
    d[strong_lookup$gene, on = setNames(c("context", "gene_id_norm"), c(context_col, ".gene_norm_tmp__")), (gene_flag_col) := 1L]
  }
  d[, .gene_norm_tmp__ := NULL]
  d
}

validate_binary_cols <- function(dt, cols, label = deparse(substitute(dt))) {
  for (col in cols) {
    vals <- dt[[col]]
    if (!all(vals %in% c(0L, 1L)) || any(is.na(vals))) {
      stop(sprintf("%s %s must be non-missing 0/1.", label, col))
    }
  }
  invisible(dt)
}

rename_eqtl_excel_cols <- function(dt) {
  d <- copy(dt)
  rename_map <- c(
    context = "expression_context",
    phenotype_id = "ENSEMBL",
    source = "eQTL_source",
    rsID = "rsid",
    GWAS = "risk_SNP_GWAS",
    DGE = "SCZD_DEG_layer_adjusted",
    DGEp05 = "SCZD_DEG_layer_adjusted_p05",
    var_coloc = "strong_coloc_variant",
    gene_coloc = "strong_coloc_gene"
  )
  old <- intersect(names(rename_map), names(d))
  if (length(old) > 0L) setnames(d, old, rename_map[old])
  d
}

rename_coloc_excel_cols <- function(dt) {
  d <- copy(dt)
  rename_map <- c(
    ds = "expression_context",
    fid = "ENSEMBL",
    cat = "hypothesis_category",
    GWAS = "risk_SNP_GWAS",
    DGE = "SCZD_DEG_layer_adjusted",
    DGEp05 = "SCZD_DEG_layer_adjusted_p05",
    eQTL = "eQTL_independent",
    eQTLnom = "eQTL_nominal"
  )
  old <- intersect(names(rename_map), names(d))
  if (length(old) > 0L) setnames(d, old, rename_map[old])
  d
}
