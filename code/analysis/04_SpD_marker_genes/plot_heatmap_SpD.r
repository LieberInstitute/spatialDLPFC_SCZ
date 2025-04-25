# Load library ----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(here)
  library(scuttle)
  library(ComplexHeatmap)
  library(viridisLite)
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
  "processed-data/rds/04_SpD_marker_genes",
  "test_enrich_PRECAST_07.rds"
))

# save a csv file for Sang Ho
write_csv(
  spd_enrich_df,
  here(
    "processed-data/rds/04_SpD_marker_genes",
    "marker_gene_enrichment_PRECAST_07.csv"
  ),
  col_names = TRUE,
)

## Load Spatial Domain Annotate DF----
# Load SpD annotation
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))



## Choose marker genes ----
# Top 20 enrichment gene for Sang Ho
sprintf("logFC_spd%02d", 1:7) |>
  map(.f = function(.var) {
    # browser()
    spd_enrich_df |>
      arrange(pick(.var)) |>
      slice_tail(n = 20) |>
      select(gene, ensembl) |>
      mutate(
        spd = str_remove(.var, "logFC_"),
        rank = 20:1
      )
  }) |>
  list_rbind() |>
  left_join(
    spd_anno_df
  ) |>
  write_csv(
    here(
      "code/analysis/04_SpD_marker_genes",
      "top_20_enrichment_SpD_maker_genes.csv"
    ),
    col_names = TRUE,
  )




# Choose top 10 most enriched genes (of logFC) for each spd
# slct_gene_df <- paste0(
#   "logFC_",
#   spd_anno_df |> arrange(anno_lab) |> pull(spd)
# ) |>
#   map(.f = function(.var) {
#     browser()
#     spd_enrich_df |>
#       arrange(pick(.var)) |>
#       slice_tail(n = 10) |>
#       select(gene, ensembl)
#   }) |>
#   list_rbind() |>
#   distinct()

# Choose top 10 t_test enriched gene for each spd
slct_gene_df <-
  paste0(
    "t_stat_",
    spd_anno_df |> arrange(anno_lab) |> pull(spd)
  ) |>
  # sprintf("t_stat_spd%02d", 1:7) |>
  map(.f = function(.var) {
    # browser()
    spd_enrich_df |>
      arrange(pick(.var)) |>
      slice_tail(n = 10) |>
      select(gene, ensembl, t_stat = .var) |>
      mutate(spd = str_remove(.var, "t_stat_"))
  }) |>
  list_rbind() |>
  distinct()

dim(slct_gene_df)


# note: quick analysis
# TODO: remove later

# spd_enrich_df |>
#   arrange(desc(logFC_spd01)) |>
#   slice_head(n = 10) |>
#   pull(gene)





# TODO: replace actual genes here
# slct_genes <- spd_enrich_df[slct_gene_ensembl, ] |>
#   with(setNames(gene, ensembl)) # Create vector of gene symbol, named by ensembl

# length(slct_genes)
anyDuplicated(slct_gene_df)
# stopifnot(anyDuplicated(slct_genes) == 0)

# Keep unique elements of slct_genes
# slct_genes <- unique(slct_genes)



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
gene_mat <- logcounts(pb_agg_spd)[slct_gene_df$ensembl, ]

# replace sample gene emsembl with gene names
rownames(gene_mat) <- slct_gene_df$gene


# TODO: annotate spd domain name and color
# TODO: add column color annotation for spatial domain
# TODO: change the spatial domain names

# TODO: Order genes

# TODO: Order columns


## Create annotation for SpD ----


# order spatial domains by cortical layers
gene_mat <- gene_mat[, spd_anno_df$spd[order(spd_anno_df$anno_lab)]]


col_annotation <- columnAnnotation(
  spd = gene_mat |> colnames(),
  col = list(
    spd = set_names(
      Polychrome::palette36.colors(7)[seq.int(7)],
      spd_anno_df$spd |> sort()
    )
  ),
  # labels = spd_anno_df$anno_lab,
  show_legend = FALSE
)

row_annotation <- rowAnnotation(
  spd = slct_gene_df$spd,
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


# Save plot
pdf(
  here(
    "plots/04_SpD_marker_genes",
    "heatmap_spd_marker_genes_enrichment_top_10_t-stat.pdf"
  ),
  height = 10
)
Heatmap(
  scaled_gene_mat,
  # TODO: choose color palette distinguiting dx-DEG and spd marker genes
  col = viridis(100), # Any sequential color palette would do.
  name = "scaled median log2 counts\n across samples",
  bottom_annotation = col_annotation, # Add bottom annotation
  left_annotation = row_annotation,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_labels = spd_anno_df$anno_lab[
    match(
      colnames(scaled_gene_mat),
      spd_anno_df$spd
    )
  ],
  cluster_rows = FALSE, # TODO: remove this
  cluster_columns = FALSE,
  row_title = "Marker Genes",
  column_title = "Spatial Domains"
) |> print()
dev.off()

# Session info ----
session_info()
