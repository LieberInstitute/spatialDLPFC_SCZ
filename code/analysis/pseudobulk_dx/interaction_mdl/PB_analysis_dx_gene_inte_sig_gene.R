# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(limma)
  # library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(pheatmap)
})

spd_rds <- list.files(
  here(
    "processed-data", "rds", "layer_spd"
  ),
  pattern = ".rds"
)

# Dx DE analysis ----
# pdf(
#   here(
#     "plots/PB_DE",
#     "test_dx_volcano_plot.pdf"
#   )
# )
.file <- "test_spe_pseudo_PRECAST_07.rds"
{
  # for (.file in spd_rds) {
  # .file <- spd_rds[1]

  .spd <- str_remove(.file, "test_spe_pseudo_") |>
    str_remove(".rds")

  sce_pseudo <- readRDS(
    here(
      "processed-data", "rds", "layer_spd",
      .file
    )
  )

  sce_pseudo$spd <- sce_pseudo[[.spd]]


  dx_mod <- model.matrix(
    ~ 0 + dx * spd + age + sex + slide_id,
    colData(sce_pseudo)
  )

  int_terms <- grep(":", dx_mod |> colnames(), value = TRUE)


  dx_sig_gene <- read_csv(
    "~/Downloads/test_PRECAST_07.csv"
  ) |> filter(fdr_scz <= 0.10)


  corfit <- duplicateCorrelation(
    logcounts(sce_pseudo)[dx_sig_gene$ensembl, ],
    design = dx_mod,
    block = colData(sce_pseudo)$sample_id
  )


  # Naive without adjusting for sample-id random effect
  fit <- lmFit(
    logcounts(sce_pseudo)[dx_sig_gene$ensembl, ],
    design = dx_mod,
    block = colData(sce_pseudo)$sample_id,
    correlation = corfit$consensus
  )
  fit <- eBayes(fit)

  Int_df <- topTable(fit, coef = c(int_terms), num = Inf)

  # cont.dif <- makeContrasts(
  #   spd01 = dxscz - dxntc,
  #   # spd02 = (dxscz + dxscz:spdspd02) - (dxntc),
  #   levels = dx_mod
  # )

  cont.mat <- rbind(
    rep(-1, 7),
    rep(1, 7),
    matrix(0, nrow = 23, ncol = 7),
    cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
  )

  colnames(cont.mat) <- sprintf("spd%02d", 1:7)

  contrast_fit <- contrasts.fit(fit, cont.mat)
  contrast_fit <- eBayes(contrast_fit)


  #
  topTable(contrast_fit) |> View()

  cont_df <- topTable(contrast_fit, coef = sprintf("spd%02d", 1:7), num = Inf)
  # colnames(cont_df) <- paste0(colnames(cont_df), "_contrast")
  cont_df <- cont_df |> rownames_to_column("gene_id")

  # saveRDS(
  #   fit,
  #   here(
  #     "processed-data/PB_dx_genes",
  #     "test_inter_PRECAST_07_sig_gene_only.rds"
  #   )
  # )

  # Omnibus test for
  # topTable(fit, coef = c("dxscz", int_terms), num = 30)

  # Omnibus test for interaction terms

  # # dx_df <- topTable(fit, coef = "dxscz", num = Inf)
  # Int_df <- topTable(fit, coef = c(int_terms), num = Inf)

  # tmp_int_df <- Int_df |> rownames_to_column("gene_id")


  # matched_df <- full_join(
  #   tmp_int_df,
  #   cont_df
  # )

  # matched_df |> select(starts_with("P.Value"))

  # int_dx_sig_gene <- dx_df |> filter(adj.P.Val <= 0.05)

  # intersect(int_dx_sig_gene$gene_name,
  # int_sig_df$gene_name
  # )

  # int_sig_df <- Int_df |>
  #   filter(adj.P.Val <= 0.15) |>
  # rownames_to_column(var = "gene_id") |>
  # left_join(
  #   rowData(sce_pseudo) |> data.frame() |> select(gene_id, gene_name)
  # )
  # rownames(int_sig_df) <- NULL
  # int_sig_df <- int_sig_df |> column_to_rownames("gene_name")

  # pheatmap(
  #   Int_df |> select(starts_with("dxscz.")) |> data.matrix(),
  #   scale = "row"
  # )

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

gene_names_hc_ordered <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)



  pdf(here("plots/PB_dx_genes",
  "dx_sig_gene_layer_specific_heatmap.pdf"),
  height = 20)
  pheatmap(
    heatmap_mat[
      gene_names_hc_ordered,
       order(colnames(heatmap_mat))],
    # scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = all_gene_mat |> transmute(`-log10P` = -1 * log10(P.Value))
  )
  dev.off()


  all_gene_mat |>
    filter(adj.P.Val <= 0.05) |>
    rownames()


  # intersect(dx_sig_gene$gene, int_sig_df$gene_name)


  # test of the deviation (FALSE)
  topTable(fit,
    coef = int_terms[6] # , n=Inf
  )

  # Test for diagnosis
  topTable(fit, coef = "dxscz", n = 30)


  # TODO: add one that adjusts effect


  #   dx_block_cor <-
  #     registration_block_cor(
  #       sce_pseudo,
  #       registration_model = dx_mod,
  #       var_sample_id = "sample_id"
  #     )

  #   dx_res <- registration_stats_enrichment(
  #     sce_pseudo,
  #     block_cor = dx_block_cor,
  #     covars = c("spd", "age", "sex", "slide_id"),
  #     var_registration = "dx",
  #     gene_ensembl = "gene_id",
  #     gene_name = "gene_name"
  #   )

  #   dir.create(
  #     here("processed-data/PB_dx_genes"),
  #     showWarnings = FALSE
  #   )

  #   ## Save DE results ----
  #   dx_res |>
  #     arrange(fdr_scz) |>
  #     write.csv(
  #       here(
  #         "processed-data/PB_dx_genes",
  #         paste0("test_", .spd, ".csv")
  #       )
  #     )

  #   ## Volcano Plot ----
  #   out_gene_df <- dx_res |>
  #     arrange(fdr_scz) |>
  #     slice_head(n = 1)
}
# dev.off()

# Session Info ----
sessioninfo::session_info()
