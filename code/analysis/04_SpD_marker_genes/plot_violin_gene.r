# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(scater)
  library(sessioninfo)
})

# Load data -----
# load psueobulked data
## SpD07 ----
.spd <- "PRECAST_07"

sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

# Select gene ----
gene_name <- "RELN"
gene_ensembl <- rowData(sce_pseudo)$gene_id[rowData(sce_pseudo)$gene_name == gene_name]

# error prevention
stopifnot(legnth(gene_ensembl) == 1)

# Plot ----
## Across all spatial domains -----

plotExpression(
  sce_pseudo,
  features = gene_ensembl,
  x = .spd,
  color = "dx",
  show_median = TRUE
) + scale_x_discrete(
  limits = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4))
  # TODO: add the layer annotation here
) +
  scale_color_manual(
    name = "Diagnosis",
    values = c(
      "ntc" = "blue", "scz" = "red"
    ),
    labels = c("NTC", "SCZ")
  )



## Between SpD01 (L6/WM) and SpD04 (WM)
# TODO: subset to SpD01 and SPD04 to run the lines above


## Oligodentrocytes ----
# from PEC data
Oligo_genes <- c("FOLH1", "CD22", "SLC5A11", "MAG", "MYRF", "ANLN", "MOG", "LDB3")
Oligo_ensembl <- rowData(sce_pseudo)$gene_id[rowData(sce_pseudo)$gene_name %in% Oligo_genes]

stopifnot(length(Oligo_genes) == length(Oligo_ensembl))

plotExpression(
  sce_pseudo,
  features = Oligo_genes,
  x = .spd,
  color = "dx",
  show_median = TRUE,
  scales = "free_y",
  swap_rownames = "gene_name"
) + scale_x_discrete(
  limits = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4))
  # TODO: add the layer annotation here
) +
  scale_color_manual(
    name = "Diagnosis",
    values = c(
      "ntc" = "blue", "scz" = "red"
    ),
    labels = c("NTC", "SCZ")
  )


## Excit L6 ----
ex_l6_genes <- c(
  # "LINC02718",
  "MCTP2",
  # "AC006299.1",
  # "DPP4",
  # "ADAMTSL1",
   "FILIP1",
    "OLFML2B", "PTPRU",
  # "AC011246.1",
  "AC024610.2"
)
# ex_l6_ensembl <- rowData(sce_pseudo)$gene_id[rowData(sce_pseudo)$gene_name %in% ex_l6_genes]

# stopifnot(length(ex_l6_genes) == length(ex_l6_ensembl))

plotExpression(
  sce_pseudo,
  features = ex_l6_genes,
  x = .spd,
  color = "dx",
  show_median = TRUE,
  scales = "free_y",
  swap_rownames = "gene_name"
) + scale_x_discrete(
  limits = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4))
  # TODO: add the layer annotation here
) +
  scale_color_manual(
    name = "Diagnosis",
    values = c(
      "ntc" = "blue", "scz" = "red"
    ),
    labels = c("NTC", "SCZ")
  )



# Save plot -----


# Session Info ----
session_info()
