# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
fld_h5ad <- here("processed-data/17_gsMap/h5ad")
dir.create(fld_h5ad, showWarnings = FALSE, recursive = TRUE)

# colData kept in the h5ad `obs`; `spd_label` is the gsMap --annotation column
keep_col <- c("key", "sample_id", "brnum", "dx", "age", "sex", "spd_label")


# Load data ----
## Finalized spe (63 donors, kept/QC-passed spots only) ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

# error prevention
stopifnot(all(keep_col %in% names(colData(spe))))
stopifnot(!anyNA(spe$spd_label))
stopifnot(!anyDuplicated(spe$key))


# Gene names ----
# gsMap maps genes to SNPs by matching gene *symbols* against its GTF
# (gencode v46lift37), so rows must be gene symbols rather than Ensembl IDs.
# For the few symbols mapping to multiple Ensembl IDs, keep the first one
# (make.unique() suffixes like "GENE.1" would not match the GTF anyway).
gene_name <- rowData(spe)$gene_name
keep_gene <- !is.na(gene_name) & !duplicated(gene_name)
cat(sprintf(
  "Dropping %d genes with missing/duplicated symbols; keeping %d genes\n",
  sum(!keep_gene), sum(keep_gene)
))
spe <- spe[keep_gene, ]
rownames(spe) <- rowData(spe)$gene_name


# Export one h5ad per sample ----
## Each Visium capture area is one gsMap "slice". gsMap needs:
##   - layers['count']: raw UMI counts (--data_layer count)
##   - obsm['spatial']: spot coordinates
##   - obs['spd_label']: annotation used by the GNN and Cauchy combination
##   - obs_names = spe$key, so per-spot results can be joined back to spe
## zellkonverter puts X_name in X and all other assays in `layers`.
export_sample_h5ad <- function(spe, .smp, fld_out) {
  sub_spe <- spe[, spe$sample_id == .smp]

  col_df <- colData(sub_spe)[, keep_col] |>
    as.data.frame() |>
    mutate(across(where(is.factor), as.character))
  rownames(col_df) <- sub_spe$key

  sce <- SingleCellExperiment(
    assays = list(
      logcounts = as(logcounts(sub_spe), "dgCMatrix"),
      count = as(counts(sub_spe), "dgCMatrix")
    ),
    colData = DataFrame(col_df),
    reducedDims = list(spatial = spatialCoords(sub_spe))
  )
  colnames(sce) <- sub_spe$key

  out_file <- file.path(fld_out, paste0(.smp, ".h5ad"))
  writeH5AD(sce, file = out_file, X_name = "logcounts")
  out_file
}

sample_ids <- sort(unique(spe$sample_id))

# error prevention
stopifnot(length(sample_ids) == 63)

h5ad_files <- sample_ids |>
  map_chr(.f = function(.smp) {
    cat("Exporting", .smp, "\n")
    export_sample_h5ad(spe, .smp, fld_h5ad)
  })

# error prevention
stopifnot(all(file.exists(h5ad_files)))


# Write sample lists used by the SLURM array jobs ----
## Line N of sample_list.txt is processed by SLURM_ARRAY_TASK_ID = N.
writeLines(sample_ids, file.path(fld_h5ad, "sample_list.txt"))
writeLines(h5ad_files, file.path(fld_h5ad, "h5ad_list.txt"))

## Per-diagnosis sample lists, for cross-sample Cauchy combination (09) ----
sample_dx <- colData(spe) |>
  as.data.frame() |>
  distinct(sample_id, brnum, dx) |>
  arrange(sample_id)

# error prevention: one diagnosis per sample
stopifnot(!anyDuplicated(sample_dx$sample_id))

writeLines(
  sample_dx |> filter(dx == "ntc") |> pull(sample_id),
  file.path(fld_h5ad, "sample_list_ntc.txt")
)
writeLines(
  sample_dx |> filter(dx == "scz") |> pull(sample_id),
  file.path(fld_h5ad, "sample_list_scz.txt")
)
write_csv(sample_dx, file.path(fld_h5ad, "sample_dx.csv"))

print("Finished exporting per-sample h5ad files for gsMap")


# Session info ----
sessioninfo::session_info()
