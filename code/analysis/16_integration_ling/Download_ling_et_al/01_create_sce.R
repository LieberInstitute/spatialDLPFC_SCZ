# Build a SingleCellExperiment from the Ling et al. (2024) donor-assigned
# count matrices (level 2: `<lib>...gmg.donors.h5`, see README.md),
# downsampled to 3 batches from 2021 = 59 donors (30 SCZ + 29 NTC).
#   - counts             : genes x nuclei, exonic + intronic UMI counts,
#                          HDF5-backed (TENx sparse format) to fit in memory
#   - colData(sce)       : per-nucleus metadata (BA46.jointMetadata.txt.gz),
#                          then donor metadata (SZvillage_donorMetadata.txt)
# Only nuclei with a joint-metadata row are kept, since the other nuclei
# carry no donor or cell-type label.
#
# Output (keep the two files in the same folder):
#   processed-data/rds/16_integration_ling/sce_ling2024_3batches/
#     - counts.h5 : TENx-format counts
#     - sce.rds   : SCE with colData/rowData but no assay; load with
#                   `load_ling_sce()` below (also defined in 02_*.R)

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(rhdf5)
  library(Matrix)
  library(HDF5Array)
  library(SingleCellExperiment)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
## Three 2021 batches within 104 days that share no donor with each other,
## so every donor is one sample in one batch: 05-21-2021 (10 SCZ + 10 NTC),
## 05-28-2021 (10 + 9) and 09-02-2021 (10 + 10). Any set of 3 batches with
## exactly 10 + 10 each and no shared donor must mix 2019 and 2021 (>= 611
## days apart), because the 2019 batches re-use donors (see README.md).
## Note 05-21-2021 used a stricter nuclei cutoff (intronic 0.300 vs 0.200).
batches_keep <- c("05-21-2021", "05-28-2021", "09-02-2021")

fld_raw <- here("raw-data/ling_et_al")
fld_out <- here("processed-data/rds/16_integration_ling/sce_ling2024_3batches")
dir.create(fld_out, showWarnings = FALSE, recursive = TRUE)
path_h5_out <- file.path(fld_out, "counts.h5")
path_sce_out <- file.path(fld_out, "sce.rds")

path_donor_meta <- file.path(fld_raw, "SZvillage_donorMetadata.txt")
path_cell_meta <- file.path(fld_raw, "BA46.jointMetadata.txt.gz")

## Per-barcode metadata (open access, ~73 MB), carries the donor assignment
## of each nucleus; the h5 files only hold barcodes.
if (!file.exists(path_cell_meta)) {
  options(timeout = 3600)
  download.file(
    paste0(
      "https://data.nemoarchive.org/other/grant/broad/mccarroll/",
      "transcriptome/sncell/10x_v3.1/human/processed/other/",
      "BA46.jointMetadata.txt.gz"
    ),
    path_cell_meta,
    mode = "wb"
  )
}


# Per-nucleus metadata ----
## Text columns (e.g. `subtype`) start with NA, so guess types from all rows
cell_meta <- read_tsv(
  path_cell_meta,
  col_types = cols(.default = col_guess(), PREFIX = col_character(),
                   CELL_BARCODE = col_character(), DONOR = col_character()),
  guess_max = Inf,
  progress = FALSE
) |>
  mutate(batch = str_extract(PREFIX, "\\d{2}-\\d{2}-\\d{4}")) |>
  filter(batch %in% batches_keep) |>
  mutate(key = paste(PREFIX, CELL_BARCODE, sep = "_"))

stopifnot(!anyDuplicated(cell_meta$key))


# Donor metadata ----
## Column names contain spaces ("Batch 1"); make them syntactic ("Batch_1")
donor_meta <- read_tsv(
  path_donor_meta,
  col_types = cols(Donor = col_character()),
  progress = FALSE
) |>
  rename_with(\(x) gsub(" ", "_", x))

stopifnot(!anyDuplicated(donor_meta$Donor))
stopifnot(all(cell_meta$DONOR %in% donor_meta$Donor))


# Check the donor design ----
donor_batch <- cell_meta |>
  distinct(batch, DONOR) |>
  left_join(donor_meta, by = c("DONOR" = "Donor"))

## Each donor in exactly one kept batch, both dx groups in every batch
stopifnot(!anyDuplicated(donor_batch$DONOR))
dx_tab <- table(donor_batch$batch, donor_batch$Schizophrenia)
dx_tab
stopifnot(all(dx_tab > 0))
message(nrow(donor_batch), " donors in ", length(batches_keep), " batches")


# Locate h5 files ----
## Skip `.part` files left by incomplete downloads
h5_files <- list.files(
  fld_raw,
  pattern = "\\.gmg\\.donors\\.h5$",
  full.names = TRUE
)
## Library name = file prefix without "_MyelinRemoved" (05-07-2019 batch),
## to match `PREFIX` in the joint metadata
names(h5_files) <- basename(h5_files) |>
  sub(pattern = "\\.intronic_.*$", replacement = "") |>
  sub(pattern = "_MyelinRemoved$", replacement = "")
h5_files <- h5_files[names(h5_files) %in% cell_meta$PREFIX]

## All 8 reactions per batch should be downloaded
stopifnot(setequal(names(h5_files), unique(cell_meta$PREFIX)))
message("Loading ", length(h5_files), " libraries from batches ",
        paste(batches_keep, collapse = ", "))


# Genes and nuclei to keep ----
## Each file only lists genes detected in that library, so a gene missing
## from a file has zero counts there. Use the union of genes across files.
## Reading only the names here lets us size the HDF5 file before writing.
lib_genes <- map(h5_files, \(f) {
  as.character(h5read(f, "matrix/features/name"))
})
lib_keys <- imap(h5_files, \(f, lib) {
  paste(lib, h5read(f, "matrix/barcodes"), sep = "_")
})
h5closeAll()

stopifnot(!any(map_lgl(lib_genes, \(g) anyDuplicated(g) > 0)))
gene_names <- sort(Reduce(union, lib_genes))
message("Genes: ", length(gene_names), " in union, ",
        length(Reduce(intersect, lib_genes)), " shared by all")

## Nuclei without a joint-metadata row have no donor, so drop them
lib_keep <- map(lib_keys, \(k) k %in% cell_meta$key)
n_drop <- sum(map_int(lib_keep, \(x) sum(!x)))
message(n_drop, " of ", sum(lengths(lib_keep)),
        " nuclei have no joint metadata entry and are dropped")
cell_keys <- unlist(map2(lib_keys, lib_keep, \(k, x) k[x]), use.names = FALSE)


# Write counts to HDF5, one library at a time ----
## Three batches (24 libraries) hold ~1.5 billion non-zero counts, near the
## dgCMatrix limit (2^31 - 1) and too large to cbind() in memory.
read_ling_h5 <- function(path, lib) {
  grp <- h5read(path, "matrix")
  h5closeAll()
  sparseMatrix(
    i = as.integer(grp$indices) + 1L,
    p = as.integer(grp$indptr),
    x = as.numeric(grp$data),
    dims = as.integer(grp$shape),
    dimnames = list(
      as.character(grp$features$name),
      paste(lib, grp$barcodes, sep = "_")
    ),
    repr = "C"
  )
}

expand_rows <- function(m, genes) {
  new("dgCMatrix",
      i = match(rownames(m), genes)[m@i + 1L] - 1L,
      p = m@p,
      x = m@x,
      Dim = c(length(genes), ncol(m)),
      Dimnames = list(genes, colnames(m))) |>
    ## row indices must be sorted within each column
    as("TsparseMatrix") |>
    as("CsparseMatrix")
}

if (file.exists(path_h5_out)) file.remove(path_h5_out)
sink <- TENxRealizationSink(
  dim = c(length(gene_names), length(cell_keys)),
  dimnames = list(gene_names, cell_keys),
  type = "integer",
  filepath = path_h5_out,
  group = "counts"
)

## The TENx sink must be filled in column order
col_start <- 1L
for (lib in names(h5_files)) {
  m <- read_ling_h5(h5_files[[lib]], lib)[, lib_keep[[lib]], drop = FALSE] |>
    expand_rows(genes = gene_names)
  stopifnot(identical(colnames(m), lib_keys[[lib]][lib_keep[[lib]]]))

  vp <- ArrayViewport(
    dim(sink),
    IRanges(c(1L, col_start), width = dim(m))
  )
  write_block(sink, vp, as(m, "SVT_SparseMatrix"))
  message(lib, ": ", ncol(m), " nuclei, ", length(m@x), " non-zero counts")

  col_start <- col_start + ncol(m)
  rm(m)
  gc(verbose = FALSE)
}
close(sink)
stopifnot(col_start - 1L == length(cell_keys))


# colData ----
col_df <- tibble(key = cell_keys) |>
  left_join(
    cell_meta |>
      rename(library = PREFIX, barcode = CELL_BARCODE),
    by = "key"
  ) |>
  left_join(donor_meta, by = c("DONOR" = "Donor")) |>
  ## Convenience columns for downstream DE
  mutate(
    dx = factor(
      if_else(Schizophrenia == "Affected", "SCZ", "NTC"),
      levels = c("NTC", "SCZ")
    ),
    batch = factor(batch, levels = batches_keep)
  ) |>
  relocate(library, barcode, batch, DONOR, dx, cellclass, subclass, subtype,
           .after = key)
stopifnot(identical(col_df$key, cell_keys))


# Build sce ----
## The counts are not stored in the rds (an HDF5-backed assay saved with
## saveRDS() keeps an absolute path); `load_ling_sce()` re-attaches them.
sce <- SingleCellExperiment(
  colData = DataFrame(select(col_df, -key), row.names = col_df$key),
  rowData = DataFrame(gene_name = gene_names, row.names = gene_names)
)
metadata(sce)$source <- list(
  study = "Ling et al. (2024) Nature 627:604-611",
  doi = "10.1038/s41586-024-07109-5",
  batches = batches_keep,
  h5_files = basename(h5_files),
  cell_metadata = basename(path_cell_meta),
  donor_metadata = basename(path_donor_meta),
  counts_h5 = basename(path_h5_out)
)
saveRDS(sce, path_sce_out)


# Check: re-load with counts ----
load_ling_sce <- function(dir) {
  sce <- readRDS(file.path(dir, "sce.rds"))
  counts <- TENxMatrix(file.path(dir, "counts.h5"), group = "counts")
  stopifnot(identical(dimnames(counts), dimnames(sce)))
  assay(sce, "counts") <- counts
  sce
}

sce <- load_ling_sce(fld_out)
sce
table(sce$batch, sce$dx)
table(sce$cellclass, sce$dx, useNA = "ifany")
colData(sce) |>
  as_tibble() |>
  distinct(batch, DONOR, dx, Sex, Age) |>
  count(batch, dx, Sex)


# Session info ----
session_info()
