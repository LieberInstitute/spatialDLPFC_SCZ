all_gene_mat <- cont_df |>
    left_join(
      rowData(sce_pseudo) |> data.frame() |> select(gene_id, gene_name)
    ) |>
    column_to_rownames("gene_name")

  # all_gene_mat <- cont_df |>
  #   rownames_to_column(var = "gene_id") |>
  #   left_join(
  #     rowData(sce_pseudo) |> data.frame() |> select(gene_id, gene_name)
  #   )
  # rownames(all_gene_mat) <- NULL
  # all_gene_mat <- all_gene_mat |> column_to_rownames("gene_name")

  # all_gene_anno_row <-

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
  pheatmap(
    heatmap_mat[, order(colnames(heatmap_mat))],
    # scale = "row",
    cluster_cols = FALSE,
    annotation_row = all_gene_mat |> transmute(`-log10P` = -1 * log10(P.Value))
  )
  dev.off()