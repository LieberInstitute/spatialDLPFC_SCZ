# Load library ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(here)
  library(scuttle)
  library(ComplexHeatmap)
  library(sessioninfo)
})

# Load data ----
## Load SpD enriched genes ----
spe_pb <- readRDS(
  here(
    # TODO: organize the rds paths
    "processed-data/rds/layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

# error prevention
stopifnot("PRECAST_07" %in% names(colData(spe_pb)))

# NOTE: code that may not that usefuls
# spe_pb$PRECAST_07 <- factor(spe_pb$PRECAST_07)
# stopifnot(is.factor(spe_pb$PRECAST_07))

## Load pseudobulked data ----
# TODO: Needs different way of showing results
spd_enrich_df <- readRDS(here(
  # TODO: need organize
  "processed-data/rds/layer_enrich_test",
  "test_enrich_PRECAST_07.rds"
))

# TODO: replace actual genes here
slct_genes <- spd_enrich_df[1:70, ] |>
  with(setNames(gene, ensembl)) # Create vector of gene symbol, named by ensembl

# Plot heatmap ----

## Aggregated across samples by median ----

# for each gene, calculate the median across all samples
pb_agg_spd <- spe_pb |>
  as("SingleCellExperiment") |>
  aggregateAcrossCells(
    ids = spe_pb$PRECAST_07,
    statistics = "median",
    use.assay.type = "logcounts"
  )

# error prevention
stopifnot("logcounts" %in% assayNames(pb_agg_spd))


## Organize data for heamtap ----
gene_mat <- logcounts(pb_agg_spd)[names(slct_genes), ]

# replace sample gene emsembl with gene names
rownames(gene_mat) <- slct_genes[rownames(gene_mat)]


# TODO: annotate spd domain name and color
# TODO: add column color annotation for spatial domain
# TODO: change the spatial domain names

# TODO: Order genes

# TODO: Order columns


## Create annotation for SpD ----
# Load SpD annotation
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

# order spatial domains by cortical layers
gene_mat <- gene_mat[, spd_anno_df$spd[order(spd_anno_df$anno_lab)]]




col_annotation <- columnAnnotation(
  spd = spd_anno_df$spd,
  col = list(
    spd = set_names(
      Polychrome::palette36.colors(7)[seq.int(7)],
      spd_anno_df$spd |> sort()
    )
  ),
  # labels = spd_anno_df$anno_lab,
  show_legend = FALSE
)


# colnames(gene_mat) <- spd_anno_df$anno_lab[
#   match(colnames(gene_mat), spd_anno_df$spd)
# ]


scaled_gene_mat <- gene_mat |>
  apply(
    MARGIN = 1, # Row-wise
    FUN = scale,
    center = TRUE,
    scale = TRUE
  ) |>
  t()

# Ensure the scaled matrix retains the column names
colnames(scaled_gene_mat) <- colnames(gene_mat)


# Create complex heatmap with annotations
# NOTE: maybe need to swap the row and columns


library(viridisLite)
pdf(
  here("test.pdf")
)
Heatmap(
  scaled_gene_mat,
  # TODO: choose color palette distinguiting dx-DEG and spd marker genes
  col = viridis(100), # Any sequential color palette would do.
  name = "scaled median log2 counts\n across samples",
  bottom_annotation = col_annotation, # Add bottom annotation
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_labels = spd_anno_df$anno_lab,
  column_order = order(spd_anno_df$anno_lab), # NOTE: Boyi think the column_order is extremely flowed.
  cluster_rows = TRUE, # TODO: remove this
  cluster_columns = FALSE,
  row_title = "Marker Genes",
  column_title = "Spatial Domains"
) |> print()
dev.off()
# ## Create complex heatmap ----
# heatmap <- pheatmap::pheatmap(
#   gene_mat,
#   scale = "row", # Scale across spd
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   # fontsize_row = 10,
#   # fontsize_col = 10,
#   legend = TRUE,
#   main = "scaled median log2 counts\n across samples"
# )

## Annotate SpD row ----


## Save plot ----



# Session info ----
session_info()
