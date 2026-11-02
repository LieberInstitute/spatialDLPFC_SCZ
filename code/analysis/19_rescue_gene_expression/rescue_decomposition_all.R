# rescue_decomposition_all.R
#
# Full spot-level decomposition for all samples: for every sample, a genes x spots
# matrix for each of the 29 cell types (alpha_kj * x_k) and one for the unattributed
# expression (r_j), i.e. a sample x gene x spot x (cell type + unattributed) object.
#
# Each rescue_<sample>.rds is a list with
#   genes       gene_id (Ensembl; rownames of all matrices) and gene_name
#   spots       spot_id, barcode, array_row, array_col, domain (fnl_spd), nUMI, unattributed_ratio
#   X           genes x K reference profiles, CPM (each column sums to 1e6); K = 29 (cellType_hc)
#   alpha.hat   K x spots, estimated cell-type weights
#   R.count     genes x spots, unattributed expression, count scale
#   Y.count     genes x spots, observed counts for the same genes and spots
#   contrib     (K + 1) x spots, read contributions (= the proportions in spot_level_results.csv)
# Model, for spot j on the CPM scale:  y_j ~ sum_k alpha_kj * x_k + r_j
#
# Storage: the matrices are genes x spots x 8 bytes each (e.g. 18,000 x 4,400 ~ 600 MB),
# so a sample takes roughly 18 GB depending on its number of spots, and all 63 samples ~1.1 TB.
# out_dir needs that much free space.
# The matrices are written one per file, so a subset of samples or cell types can be run instead.
#   <out_dir>/<sample>/<cell type>.rds    genes x spots
#   <out_dir>/<sample>/Unattributed.rds   genes x spots
#   <out_dir>/<sample>/genes.csv, spots.csv, cell_types.csv (cell type -> file name)

library(Matrix)

in_dir  <- "./rescue_out"     # folder with the rescue_<sample>.rds files
out_dir <- "decomposition"
scale   <- "count"            # "count" (same scale as Y.count) or "cpm" (RESCUE's L.hat / R.hat scale)

files <- list.files(in_dir, pattern = "^rescue_.*\\.rds$", full.names = TRUE)
fname <- function(x) gsub("[^[:alnum:]_.-]", "_", x)

for (f in files[1]) {            # this requires large  per sample
  res <- readRDS(f)
  X <- res$X; alpha <- res$alpha.hat
  nUMI <- colSums(res$Y.count)
  sc <- if (scale == "count") nUMI / 1e6 else rep(1, length(nUMI))   # alpha_kj * x_k is CPM; * nUMI_j / 1e6 gives counts
  
  d <- file.path(out_dir, res$sample_id)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (k in colnames(X)) {
    A <- outer(X[, k], alpha[k, ] * sc)                                # genes x spots
    saveRDS(A, file.path(d, paste0(fname(k), ".rds")), compress = FALSE)
  }
  R <- if (scale == "count") res$R.count else res$R.count %*% Diagonal(x = 1e6 / nUMI)
  dimnames(R) <- dimnames(res$R.count)
  saveRDS(R, file.path(d, "Unattributed.rds"))
  
  write.csv(res$genes, file.path(d, "genes.csv"), row.names = FALSE)
  write.csv(res$spots, file.path(d, "spots.csv"), row.names = FALSE)
  write.csv(data.frame(cell_type = c(colnames(X), "Unattributed"), file = c(paste0(fname(colnames(X)), ".rds"), "Unattributed.rds")),
            file.path(d, "cell_types.csv"), row.names = FALSE)
  cat(res$sample_id, ":", ncol(X), "cell types + unattributed,", nrow(X), "genes x", ncol(alpha), "spots\n")
  rm(res, A, R); gc()
}

# Reading back. One sample as a genes x spots x (K + 1) array (needs ~18 GB of memory):
#   s  <- "V12D07-334_A1"
#   ct <- read.csv(file.path(out_dir, s, "cell_types.csv"))
#   arr <- simplify2array(lapply(file.path(out_dir, s, ct$file), function(f) as.matrix(readRDS(f))))
#   dimnames(arr)[[3]] <- ct$cell_type
