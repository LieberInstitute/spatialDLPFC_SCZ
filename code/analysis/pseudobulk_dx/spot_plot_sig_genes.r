suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(readxl)
  library(escheR)
  library(sessioninfo)
})

# Load Data ----
## Load SpatialClustering SPE Obejct ----
raw_spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)


# raw_spe <- readRDS(
#   here(
#     "processed-data/rds/02_visium_qc",
#     "qc_spe_wo_spg_N63.rds"
#   )
# )



# Subsetting two samples -
sample_names <- c(
  "V13M06-342_D1",
  "V13M06-343_D1"
)

spe <- raw_spe[, spe$sample_id %in% sample_names]

## Load Excel File ----
gene_names <- readxl::read_excel(
  here(
    "code/analysis/visium_spatial_clustering",
    "Escher_genelist.xlsx"
  ),
  col_names = FALSE
) |> unlist()
names(gene_names) <- NULL


# Run log-normaliztion for the data -----
if (!"logcounts" %in% assayNames(spe)) {
  spe <- scater::logNormCounts(
    spe,
    size.factors = sizeFactors(spe),
    transform = "log"
  )
}



## FINDING the ensemble of these datasets

gene_ensembl <- rowData(spe)$gene_id[match(gene_names, rowData(spe)$gene_name)]

## Copy the gene to colData
logcounts_mat <- logcounts(spe)[gene_ensembl, ]

rownames(logcounts_mat) <- gene_names


stopifnot(
  identical(
    colnames(logcounts_mat),
    rownames(colData(spe))
  )
)

colData(spe) <- cbind(
  colData(spe),
  logcounts_mat |> t() |> data.matrix()
)


NTC_spe <- spe[, spe$sample_id == "V13M06-342_D1"]
SCZ_spe <- spe[, spe$sample_id == "V13M06-343_D1"]

## Concatenate matrix to colData
# TODO: check if the column name is the gene name instead of ensembl_id.


# Plot indiviudal genes ----
for (gene in gene_names) {
  # browser()
  ## escheR of Sample_1
  NTC_plot <- NTC_spe |>
    make_escheR() |>
    add_ground(var = "PRECAST_07") |>
    add_fill(var = gene) +
    scale_fill_gradient(low = "white", high = "black")

  ## escheR of Sample_2
  SCZ_plot <- SCZ_spe |>
    make_escheR() |>
    add_ground(var = "PRECAST_07") |>
    add_fill(var = gene) +
    scale_fill_gradient(low = "white", high = "black")

  spot_plot <- ggpubr::ggarrange(NTC_plot, SCZ_plot,
    common.legend = TRUE
  )

  # print(spot_plot)

  pdf(
    here(
      "plots/PB_dx_genes/spot_plot",
      paste0("spot_plot_", gene, ".pdf")
    ),
    width = 8, height = 6
  )
  print(spot_plot)
  dev.off()
}


# Session Info----
session_info()
