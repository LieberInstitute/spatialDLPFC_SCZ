library(here)
library(SingleCellExperiment)
library(DropletUtils)
library(tidyverse)
library(sessioninfo)

# NOTE: to load the snRNA data in, we need to set the word directory
setwd("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/")
sce <- readRDS("se.rds")

stopifnot(
  here() == "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100"
)



# Down sampling the size ----
table(sce$cellType_broad_hc)
# Remove ambiguous cell type ----
sce <- sce[, !sce$cellType_broad_hc %in% c("Ambiguous")]

# Subset samples
table(sce$BrNum)
unique(sce$BrNum) |> length()
unique(sce$SAMPLE_ID) |> length()
## Overlapping samples
sub_samples <- readxl::read_xlsx(
  here::here("code/xenium_panel_design/Xenium_DonorList_Edit.xlsx"),
  col_names = FALSE
)[, 1:2] |> unlist()

overlap_samples <- intersect(unique(sce$BrNum), sub_samples)
# [1] "Br6432" "Br8667"
sub_sce <- sce[, sce$BrNum %in% overlap_samples]
ncol(sub_sce)
table(sub_sce$cellType_broad_hc)

set.seed(20240827)
subset_ind <- c(
  # Downsample from cell types whose counts are more than 5000
  sample(which(sce$cellType_broad_hc == "Oligo"), size = 5000),
  sample(which(sce$cellType_broad_hc == "Excit"), size = 5000),
  sample(which(sce$cellType_broad_hc == "Inhib"), size = 5000),
  # all cells from cell types whose counts are less than 5000
  which(sce$cellType_broad_hc %in% c("Astro", "EndoMural", "Micro", "OPC")),
  # ll cells from over lapping samples
  which(sce$BrNum %in% overlap_samples)
) |>
  unique() |>
  sort()

length(subset_ind)

sce <- sce[, subset_ind]
ncol(sce)
# [1] 29854



# Calculate

# Create ref file ----
dir.create(
  here("processed-data/Xenium_probe_design"),
  showWarnings = FALSE
)


if (file.exists(here("processed-data/Xenium_probe_design/PEC_snRNA_ref"))) {
  unlink(
    here("processed-data/Xenium_probe_design/PEC_snRNA_ref"),
    recursive = TRUE
  )
}

stopifnot(!file.exists(here("processed-data/Xenium_probe_design/PEC_snRNA_ref")))

# Create file ----
## count matrix file ----
write10xCounts(
  here("processed-data/Xenium_probe_design/PEC_snRNA_ref"),
  x = as(counts(sce), "sparseMatrix"),
  gene.id = rowData(sce)$gene_id,
  gene.symbol = rownames(sce),
  barcodes = colnames(sce),
  type = "sparse",
  version = "3"
)

## Create cell annotation file ----
write.table(
  colData(sce) |>
    data.frame() |>
    transmute(
      barcode = Barcode,
      annotation = cellType_broad_hc
    ),
  file = here(
    "processed-data/Xenium_probe_design/PEC_snRNA_ref",
    "annotations.csv"
  ),
  quote = FALSE,
  sep = ",", row.names = FALSE
)


# Validation ----
list.files(here("processed-data/Xenium_probe_design/PEC_snRNA_ref"))

## Check size ----
# NOTE: should be smaller than 500 MB
system("du -sh /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Xenium_probe_design/PEC_snRNA_ref")


# Create zip file ---
utils::zip(
  here(
    "processed-data/Xenium_probe_design",
    "PEC_snRNA_ref.zip"
  ),
  list.files(
    here("processed-data/Xenium_probe_design/PEC_snRNA_ref"),
    full.names = TRUE
  ),
  zip = "zip"
)
