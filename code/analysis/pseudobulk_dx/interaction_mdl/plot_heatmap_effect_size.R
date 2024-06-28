# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
})

# Create Contrast Matrix ----
  cont.mat <- rbind(
    rep(-1, 7),
    rep(1, 7),
    matrix(0, nrow = 23, ncol = 7),
    cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
  )

  colnames(cont.mat) <- sprintf("spd%02d", 1:7)

# Calculate the contrast

fit <- readRDS(
  here(
      "processed-data/PB_dx_genes",
      "test_inter_PRECAST_07_20240627.rds"
    )
)

  contrast_fit <- contrasts.fit(fit, cont.mat)
  contrast_fit <- eBayes(contrast_fit)


cont_df <- topTable(contrast_fit, coef = sprintf("spd%02d", 1:7), num = Inf)
  # colnames(cont_df) <- paste0(colnames(cont_df), "_contrast")
  cont_df <- cont_df |> rownames_to_column("gene_id")


gene_names_hc_ordered <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)


  sce_pseudo <- readRDS(
    here(
      "processed-data", "rds", "layer_spd",
      "test_spe_pseudo_PRECAST_07.rds"
    )
  )

  all_gene_mat <- cont_df |>
    left_join(
      rowData(sce_pseudo) |> data.frame() |> select(gene_id, gene_name)
    ) |>
    filter(gene_name %in% gene_names_hc_ordered) |> 
    column_to_rownames("gene_name")

stopifnot(
  nrow(all_gene_mat) == length(gene_names_hc_ordered)
  )

  spd_anno_df <- read_csv(
    here(
      "processed-data/man_anno",
      "spd_labels_k7.csv"
    )
  ) |>
    mutate(anno_lab = paste0(label, " (", spd, ") "))


  heatmap_mat <- all_gene_mat |>
    select(starts_with("spd")) |>
    data.matrix()
  colnames(heatmap_mat) <- spd_anno_df$anno_lab[match(colnames(heatmap_mat), spd_anno_df$spd)]

 pdf(here("plots/PB_dx_genes",
  "dx_sig_gene_layer_specific_heatmap.pdf"),
  height = 20)
  Heatmap(
    heatmap_mat[
      gene_names_hc_ordered,
       order(colnames(heatmap_mat))],
    # scale = "row",
    cluster_rows = FALSE,
    cluster_columns = FALSE,

    right_annotation = rowAnnotation(`sig_value` = all_gene_mat[gene_names_hc_ordered, ] |>
     transmute(`-log10P` = -1 * log10(P.Value)) |> pull(`-log10P`)
    )
  )
  dev.off()