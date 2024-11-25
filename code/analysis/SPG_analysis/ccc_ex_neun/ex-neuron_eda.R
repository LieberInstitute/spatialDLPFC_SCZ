# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scater)
  library(tidyverse)
  library(ggbeeswarm)
  library(escheR)
  library(sessioninfo)
})

# Load spe data ----
# use V12F14-053_A1 as example
# NOTE: this test dataset already contains spatial domain & SPG calling
sub_spe <- readRDS(
  here(
    "processed-data/rds/PB_dx_spg",
    "test_small_spe_for_neighbor.rds"
  )
)

## Subset to Neun+ spots only
neun_spe <- sub_spe[, sub_spe$neun_pos == TRUE]


# Analysis ----
# NOTE: marker genes for excitatory neurons - high SLC17A7 transcripts, but with low GAD1 and GAD2 transcripts.

mk_gene_ind <- which(rowData(neun_spe)$gene_name %in% c("SLC17A7", "GAD1", "GAD2"))
stopifnot(length(mk_gene_ind) == 3)




## Viz pairwise correaltion ----
marker_gene_df <- logcounts(neun_spe)[mk_gene_ind, ] |> data.matrix()
rownames(marker_gene_df) <- rowData(neun_spe)$gene_name[mk_gene_ind]
pairs(marker_gene_df |> t())
# I think this is much more complex signal than previous spot calling based on one single marker.
# Boyi decides to use PCA, etc to see how it goes.

## TSNE viz based on three marker genes ----
set.seed(20241122)
tmp_spe <- neun_spe[mk_gene_ind, ]
tmp_spe <- runTSNE(tmp_spe, dimred = NULL, assay.type = "logcounts")
plotReducedDim(tmp_spe, dimred = "TSNE", color_by = "PRECAST_07")

## UMAP viz
set.seed(20241122)
tmp_spe <- runUMAP(tmp_spe, dimred = NULL, assay.type = "logcounts")
plotReducedDim(tmp_spe, dimred = "UMAP", colour_by = "PRECAST_07")


## Visualize cluster in spatial domain.

# Session Info ----
sessioninfo::session_info()
