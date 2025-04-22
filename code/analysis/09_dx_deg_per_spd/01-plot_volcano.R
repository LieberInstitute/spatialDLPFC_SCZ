# Plot ----
## Plot Volcano ----
# browser()
# ggplot(results_enrichment) +
#   geom_point(
#     aes(
#       x = logFC_scz,
#       y = -log10(fdr_scz))
#   ) +
#   title()
impl_gene_df <- results_enrichment |>
  filter(gene %in% c(
    "PVALB",
    "NOS1",
    "SST",
    "CHODL",
    "GRIN2A",
    "SV2A",
    "DLG4",
    "C4A",
    "C3"
  )) |>
  select(ensembl, gene, ends_with("scz"))

(ggplot(
  results_enrichment,
  aes(x = logFC_scz, y = -log10(fdr_scz))
) +
  geom_point() +
  geom_label_repel(
    data = impl_gene_df, # Add labels last to appear as the top layer
    aes(label = gene),
    force = 2,
    nudge_y = 0.1
  ) +
  labs(
    title = paste0(
      "PB-", .spd_var,
      "-",
      spe_sub$spd |> unique() |> as.character()
    )
  ) +
  theme_minimal()) |>
  print()

## Plot Gene expression as violin plot ----
remain_gene_ids <- intersect(
  c(
    "ENSG00000100362",
    "ENSG00000089250",
    "ENSG00000157005",
    "ENSG00000154645"
  ),
  rowData(sce_pseudo)$gene_id
)

ge_df <- remain_gene_ids |>
  map_dfc(~ logcounts(sce_pseudo)[.x, ])


#
# (data.frame(
#   dx = sce_pseudo$dx,
#   PVALB_log = logcounts(sce_pseudo)["ENSG00000100362",],
#   NOS1_log = logcounts(sce_pseudo)["ENSG00000089250",],
#   SST_log = logcounts(sce_pseudo)["ENSG00000157005",]#,
#   # CHODL_log = logcounts(sce_pseudo)["ENSG00000154645",]
# ) |>

if (ncol(ge_df) < 1) {
  (
    ggplot() +
      labs(
        title = paste0(
          "PB-", .spd_var,
          "-",
          spe_sub$spd |> unique() |> as.character()
        )
      )
  ) |> print()
} else {
  colnames(ge_df) <- paste0(
    rowData(sce_pseudo)[remain_gene_ids, "gene_name"],
    "_log"
  )
  ge_df <- cbind(
    ge_df,
    dx = sce_pseudo$dx
  )

  (
    ge_df |>
      pivot_longer(
        cols = ends_with("_log")
      ) |>
      ggplot() +
      geom_violin(aes(x = dx, y = value)) +
      facet_wrap(vars(name), scales = "free") +
      labs(
        title = paste0(
          "PB-", .spd_var,
          "-",
          spe_sub$spd |> unique() |> as.character()
        )
      )
  ) |> print()
}
