# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scater)
  library(tidyverse)
  library(ggbeeswarm)
  library(escheR)
  library(sessioninfo)
  library(scran)
  library(bluster)
  library(pheatmap)
  library(readxl)
})

# Load data ----
## Load SPG spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

### Call SPG + spots ----
# Call SPG spots
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)


stopifnot(is.null(spe$neun_pos) == FALSE)

## Subset to Neun+ spots only ----
neun_spe <- spe[, spe$neun_pos == TRUE]

# ## Load Maker genes excel sheet ----
mk_gene_raw <- read_xlsx(
  here(
    "code/xenium_panel_design/TableS13_marker_stats_supp_table.xlsx"
  ),
  sheet = "marker_stats_supp_table"
)

mk_gene <- mk_gene_raw |>
  dplyr::filter(cellTypeResolution == "broad") |>
  group_by(cellType.target) |>
  arrange(desc(logFC)) |>
  slice_head(n = 100) |>
  ungroup() |>
  filter(cellType.target %in% c("Excit", "Inhib"))

# mk_gene$gene |> unique()


## Before integration UMAP ----
set.seed(20241202)
neun_spe <- fixedPCA(neun_spe, subset.row = mk_gene$gene |> unique())
neun_spe <- runUMAP(neun_spe, dimred = "PCA")

plotReducedDim(neun_spe,
  dimred = "UMAP", colour_by = "sample_id",
  scattermore = TRUE
)

# NOTE: Samples looks farely integratedk

## K-means (k=2) clustering ----
clust_kmeans <- clusterCells(neun_spe,
  use.dimred = "PCA",
  BLUSPARAM = KmeansParam(centers = 2)
)
colLabels(neun_spe) <- clust_kmeans

plotReducedDim(neun_spe,
  dimred = "UMAP", colour_by = "label",
  scattermore = TRUE
)

# Heatmap for marker genes
tmp_ensmbl <- mk_gene |>
  filter(cellType.target %in% c("Excit", "Inhib")) |>
  pull(gene)

# logcounts(neun_spe)[tmp_ensmbl, ] |>
#   pheatmap(
#     scale = "row",
#     cluster_rows = TRUE,
#     cluster_cols = TRUE,
#     show_rownames = TRUE,
#     show_colnames = FALSE,
#     annotation_col = as.data.frame(colData(neun_spe)$label |> factor())
#   )

## Heatmap of marker genes of pseudobulked data
pseudobulk_sce <- scuttle::aggregateAcrossCells(
  neun_spe,
  ids = interaction(
    neun_spe$sample_id, neun_spe$label
  ),
  statistics = "sum",
  use.dimred = FALSE
)

pseudobulk_sce <- logNormCounts(
  pseudobulk_sce,
  size.factors = NULL
)

anno_col_df <- data.frame(
  labels = set_names(
    pseudobulk_sce$label |> as.character(),
    nm = pseudobulk_sce$ids |> as.character()
  )
)

anno_row_df <- mk_gene |>
  select(gene, cellType.target) |>
  column_to_rownames("gene")

pheatmap(
  mat = logcounts(pseudobulk_sce)[tmp_ensmbl, ],
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = anno_col_df,
  annotation_row = anno_row_df
)

## Marker Gene Distribution ----
# SLC17A7 - enriched in Ex cells
library(scater)
plotExpression(pseudobulk_sce, features = c("ENSG00000104888"), x = "label", color_by = "label", point_size = 0.1)

# GAD1 - depleted in Ex cells
plotExpression(pseudobulk_sce, features = c("ENSG00000128683"), x = "label", color_by = "label", point_size = 0.1)

# GAD2 - depleted in Ex cells
plotExpression(pseudobulk_sce, features = c("ENSG00000136750"), x = "label", color_by = "label", point_size = 0.1)

plotExpression(pseudobulk_sce, features = c("ENSG00000136750"), x = "label", color_by = "label", point_size = 0.1)

## Donor level ---
neun_spe$sample_label_combo <- interaction(neun_spe$sample_id, neun_spe$label)
plotExpression(neun_spe, features = c("ENSG00000104888"), x = "sample_label_combo", color_by = "label", point_size = 0.1, scattermore = TRUE)

# Conclusion: Label 1 is Excitatry cell

neun_spe$ex_cell <- neun_spe$label == 1

saveRDS(
  neun_spe,
  file = here("code/analysis/SPG_analysis/ccc_ex_neun", "ex_neun_spe.rds")
)

# Session Info ----
sessioninfo::session_info()