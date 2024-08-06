# Deprecated 
# NOTE: decide to select cell-type marker genes heuristically.

# library(SingleCellExperiment)

# # library(spatialLIBD)
# library(readxl)
# library(here)
# library(tidyverse)
# library(scran)
# library(bluster)
# library(sessioninfo)


# # Load Data ----
# ## Load Batch-corrected data with finalized cell types ----

# load(
#   file.path(
#     "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq",
#     "processed-data/sce/",
#     "sce_DLPFC.Rdata"
#   )
# )

# ls()
# # [1] "sce"

# sce$cellType_broad_hc |> table()
# # Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib Ambiguous
# #  3979      2157      1601     10894      1940     24809     11067     21157
# sce$cellType_broad_k |> table()
# #  Astro  EndoMural MicroOligo      Oligo        OPC      Excit      Inhib
# #   3557       1330       5541      33716       1791      21233      10413
# #   drop
# #     23

# ## Load Maker genes excel sheet ----
# mk_gene_raw <- read_xlsx(
#   here(
#     "code/xenium_panel_design/TableS13_marker_stats_supp_table.xlsx"
#   ),
#   sheet = "marker_stats_supp_table"
# )

# mk_gene <- mk_gene_raw |>
#   dplyr::filter(cellTypeResolution == "broad") |>
#   group_by(cellType.target) |>
#   arrange(desc(ratio)) |>
#   slice_head(n = 5)

# mk_gene_min_OPC <- mk_gene |>
#   dplyr::filter(cellType.target != "OPC")

# # Broad cell types
# mk_gene$cellType.target |> table()
# # Astro EndoMural     Excit     Inhib     Micro     Oligo       OPC
# #     5         5         5         5         5         5         5

# # Num of unique genes
# mk_gene$gene |>
#   unique() |>
#   length()


# # Define Maker genes & experiment data----
# exp_sce <- sce[mk_gene_min_OPC$symbol, ]

# # Clustering ----
# ## K-means ----
# clust.kmeans2 <- clusterCells(
#   exp_sce,
#   assay.type = "binomial_deviance_residuals",
#   BLUSPARAM = KmeansParam(centers = 7)
# )



# # Examination -----



# # Session Info ----
# session_info()
