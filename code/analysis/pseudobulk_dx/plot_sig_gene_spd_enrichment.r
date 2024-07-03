# Load packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(ComplexHeatmap)
  library(viridis)
  library(sessioninfo)
})

# Load Dx_DEG data ---
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
)

sig_gene <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    # "Test_90DEGs.xlsx"
    "Test_67DEGs.xlsx"
  ),
  col_names = FALSE
)[[1]]

sig_gene_df <- gene_df |>
  filter(gene %in% sig_gene)

n_gene <- length(sig_gene)

ann_df <- sig_gene_df |>
  column_to_rownames(var = "gene") |>
  transmute(
    SCZ_reg = factor(
      logFC_scz > 0,
      levels = c(TRUE, FALSE),
      labels = c("Up", "Down")
    )
  )

# Format SpD to readable form
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

# Load PB Data ----
## All data ---
spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

## NTC only ---

ntc_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07_ntc.rds"
  )
)

## SCZ only ---

scz_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07_scz.rds"
  )
)

# Process PB Data

## Function ---
create_median_from_pb_data <- function(sce) {
  col_df <- colData(sce) |>
    data.frame() |>
    select(PRECAST_07, dx, sample_id) |>
    mutate(
      PRECAST_07 = factor(PRECAST_07,
        levels = c(
          "spd07", "spd06", "spd02",
          "spd05", "spd03", "spd01", "spd04"
        )
      )
    ) |>
    arrange(
      PRECAST_07, dx, sample_id
    )

  gene_mat <- logcounts(sce)[sig_gene_df$ensembl, rownames(col_df)]

  # Need to accumulate genes
  gene_mat_long <- gene_mat |>
    data.frame() |>
    rownames_to_column(var = "gene") |>
    pivot_longer(
      cols = starts_with("V"),
      names_to = c("sample_id", "spd"),
      names_pattern = "(.*)_(spd.*)"
    )

  gene_mat_median <- gene_mat_long |>
    group_by(gene, spd) |>
    summarize(median = median(value)) |>
    ungroup() |>
    pivot_wider(
      id_cols = "gene",
      names_from = "spd",
      values_from = "median"
    ) |>
    column_to_rownames("gene") |>
    data.matrix()

  # Convert Ensembl to Gene name
  rownames(gene_mat_median) <- rowData(sce)[
    rownames(gene_mat_median),
    "gene_name"
  ]


  colnames(gene_mat_median) <- spd_anno_df$anno_lab[
    match(colnames(gene_mat_median), spd_anno_df$spd)
  ]

  # Code to scale each row
  # only necessary when using ComplexHeatmap::Heatmap
  gene_mat_median_scaled <- apply(
    gene_mat_median,
    MARGIN = 1,
    FUN = scale
  ) |> t()

  return(gene_mat_median)
}






# Heat map ----

all_data_median <- create_median_from_pb_data(spe_pb)
ntc_data_median <- create_median_from_pb_data(ntc_pb)
scz_data_median <- create_median_from_pb_data(scz_pb)



# hc <- hclust(dist(gene_mat_median_scaled))

# gene_names_hc_ordered <- rownames(gene_mat_median)[hc$order]

# saveRDS(
#   gene_names_hc_ordered,
#   here(
#     "code/analysis/pseudobulk_dx",
#     "spd_hierarchical_cluster_order.rds"
#   )
# )

gene_names_hc_ordered <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    sprintf(
      "spd_hierarchical_cluster_order_%02d_gene.rds",
      n_gene
    )
  )
)



heatmap_all <- ComplexHeatmap::pheatmap(
  mat = all_data_median[
    gene_names_hc_ordered,
    order(colnames(all_data_median))
  ],
  name = "ALL (Scaled median logCPM)",
  color = viridis(100, option = "magma"),
  scale = "row",
  column_title = "Enrichment",
  # cluster_rows = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  row_split = ann_df[gene_names_hc_ordered, ],
  annotation_row = ann_df[gene_names_hc_ordered, , drop = FALSE],
  show_row_dend = FALSE,
  # annotation_col = col_df |> select(
  #   PRECAST_07,
  #   # sample_id, # Overwhelm color pallete
  #   dx
  # ),
  cellwidth = 10,
  cellheight = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_colors = list(SCZ_reg = c("Up" = "red", "Down" = "blue"))
)

pdf(
  file = here(
    "plots/PB_dx_genes/",
    sprintf(
    "test_sig_gene_enrich_SCZ_spd_median_logCPM_%02dGene.pdf",
    n_gene
    )

  ),
  height = 20
)
plot(heatmap_all)
dev.off()

# gene_names_hc_ordered <- heatmap_all@row_names_param$labels[heatmap_all |> row_order() |> unlist()]

# saveRDS(
#   gene_names_hc_ordered,
#   here(
#     "code/analysis/pseudobulk_dx",
#     "spd_hierarchical_cluster_order.rds"
#   )
# )

# Session Info ----
session_info()
