# Astrocyte-specific gene expression from the RESCUE decomposition.
#
# For spot j, RESCUE models (CPM scale)  y_j ~ sum_k alpha_kj * x_k + r_j.
# The astrocyte part on the count scale (same scale as Y.count) is
#   astro_gj = X[g, "Astro"] * alpha.hat["Astro", j] * nUMI_j / 1e6
# This script builds that genes x spots matrix for all 63 samples and stores
# it as the `counts` assay of the finalized spe, keeping all colData.

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(Matrix)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
ct <- "Astro" # cell type column in RESCUE's X / alpha.hat

fld_rescue <- here("processed-data/RESCUE/rescue_out")
fld_out <- here("processed-data/rds/19_rescue_gene_expression/astrocyte")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)


# Load data ----
## Finalized spe with RESCUE proportions in colData ----
spe <- readRDS(
  here(
    "processed-data/rds/18_rescue_cell_type_proportion",
    "fnl_spe_kept_spots_only_w_RESCUE.rds"
  )
)
stopifnot(identical(colnames(spe), spe$key))

## Astrocyte profile and weights from each RESCUE output ----
# Only the small pieces are kept (profile, weights, nUMI, ids), not the
# genes x spots matrices.
rescue_files <- list.files(
  fld_rescue,
  pattern = "^rescue_.*\\.rds$", full.names = TRUE
)

rescue_ls <- map(rescue_files, \(f) {
  res <- readRDS(f)
  stopifnot(ct %in% colnames(res$X), ct %in% rownames(res$alpha.hat))
  list(
    sample_id = res$sample_id,
    gene_id = res$genes$gene_id,
    x = res$X[, ct],
    spot_id = res$spots$spot_id,
    alpha = res$alpha.hat[ct, ],
    nUMI = res$spots$nUMI
  )
}) |>
  (\(l) set_names(l, map_chr(l, "sample_id")))()

cat(sprintf("Read %d RESCUE outputs\n", length(rescue_ls)))


# Align genes and spots ----
## Genes: union of RESCUE genes across samples, in spe row order ----
# RESCUE uses a slightly different gene set per sample. A gene RESCUE did not
# use in a sample gets 0 astrocyte expression in that sample.
rescue_genes <- unique(unlist(map(rescue_ls, "gene_id")))
stopifnot(all(rescue_genes %in% rownames(spe)))
genes_keep <- rownames(spe)[rownames(spe) %in% rescue_genes]

n_smp_per_gene <- table(unlist(map(rescue_ls, "gene_id")))[genes_keep]

## Spots: one-to-one with spe ----
all_spots <- unlist(map(rescue_ls, "spot_id"))
stopifnot(
  !anyDuplicated(all_spots),
  length(all_spots) == ncol(spe),
  all(all_spots %in% colnames(spe))
)

# Fill columns in spe order, so spots must be contiguous by sample
smp_order <- unique(spe$sample_id)
stopifnot(
  all(smp_order %in% names(rescue_ls)),
  !is.unsorted(match(spe$sample_id, smp_order))
)


# Build sparse astrocyte expression matrix ----
# Assembled directly from dgCMatrix slots (i, p, x) to avoid holding
# per-sample copies plus a cbind result. Only spots with alpha > 0 and genes
# with x > 0 contribute nonzeros.
nnz_per_smp <- map_dbl(rescue_ls, \(s) sum(s$alpha > 0) * sum(s$x > 0))
nnz <- sum(nnz_per_smp)
cat(sprintf("Nonzero entries: %.0f\n", nnz))
# error prevention: dgCMatrix indexes nonzeros with 32-bit integers
stopifnot(nnz < .Machine$integer.max)

mat_i <- integer(nnz) # 0-based row index
mat_x <- numeric(nnz)
col_nnz <- integer(ncol(spe))

pos <- 0
for (s_id in smp_order) {
  s <- rescue_ls[[s_id]]

  # genes with nonzero profile, as sorted row indices of the output matrix
  g_row <- match(s$gene_id, genes_keep)
  g_ok <- s$x > 0
  g_ord <- order(g_row[g_ok])
  g_row <- g_row[g_ok][g_ord]
  g_x <- s$x[g_ok][g_ord]

  # spots in spe column order
  spe_col <- which(spe$sample_id == s_id)
  m <- match(colnames(spe)[spe_col], s$spot_id)
  stopifnot(!anyNA(m))
  w <- s$alpha[m] * s$nUMI[m] / 1e6 # alpha_j * nUMI_j / 1e6

  for (j in which(w > 0)) {
    idx <- pos + seq_along(g_row)
    mat_i[idx] <- g_row - 1L
    mat_x[idx] <- g_x * w[j]
    col_nnz[spe_col[j]] <- length(g_row)
    pos <- pos + length(g_row)
  }
  cat(s_id, ":", sum(w > 0), "of", length(w), "spots with astrocyte weight > 0\n")
}
stopifnot(pos == nnz)

astro_counts <- new(
  "dgCMatrix",
  i = mat_i, p = c(0L, cumsum(col_nnz)), x = mat_x,
  Dim = c(length(genes_keep), ncol(spe)),
  Dimnames = list(genes_keep, colnames(spe))
)
rm(mat_i, mat_x)
gc()


# Checks ----
## Recompute one sample directly from its RESCUE output ----
s_chk <- smp_order[1]
res <- readRDS(file.path(fld_rescue, paste0("rescue_", s_chk, ".rds")))
direct <- outer(res$X[, ct], res$alpha.hat[ct, ] * res$spots$nUMI / 1e6)
dimnames(direct) <- list(res$genes$gene_id, res$spots$spot_id)
built <- as.matrix(astro_counts[res$genes$gene_id, res$spots$spot_id])
stopifnot(isTRUE(all.equal(built, direct, check.attributes = FALSE)))
# genes RESCUE did not use in this sample are 0
stopifnot(sum(astro_counts[!genes_keep %in% res$genes$gene_id, res$spots$spot_id]) == 0)
rm(res, direct, built)

## Total astrocyte counts per spot = alpha_j * nUMI_j ----
# X columns sum to 1e6 over each sample's genes
alpha_spe <- unlist(map(rescue_ls, "alpha"))[match(colnames(spe), all_spots)]
nUMI_spe <- unlist(map(rescue_ls, "nUMI"))[match(colnames(spe), all_spots)]
stopifnot(isTRUE(all.equal(
  colSums(astro_counts), alpha_spe * nUMI_spe,
  check.attributes = FALSE
)))


# Build astrocyte-specific spe ----
# Drop the observed assays before subsetting rows to avoid copying them
assays(spe) <- list()
spe <- spe[genes_keep, ]
assay(spe, "counts", withDimnames = FALSE) <- astro_counts
rm(astro_counts)

rowData(spe)$RESCUE_n_samples <- as.integer(n_smp_per_gene)
spe$RESCUE_alpha_Astro <- alpha_spe
spe$astro_sum_umi <- alpha_spe * nUMI_spe

metadata(spe)$RESCUE_cell_type_GE <- list(
  cell_type = ct,
  assay = "counts: X[, ct] * alpha.hat[ct, ] * nUMI / 1e6 (expected, non-integer counts)",
  source = fld_rescue
)

print(spe)


# Save ----
saveRDS(
  spe,
  file.path(fld_out, "fnl_spe_kept_spots_only_astrocyte_specific_GE.rds")
)

print("Finished building astrocyte-specific gene expression spe")


# Session info ----
sessioninfo::session_info()
