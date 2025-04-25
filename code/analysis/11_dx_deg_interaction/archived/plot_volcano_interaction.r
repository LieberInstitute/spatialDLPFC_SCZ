# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
  library(ggrepel)
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


# tmp_genes <- topTable(contrast_fit, coef = sprintf("spd%02d", 1:7), num =10) |> rownames()

# topTable(contrast_fit, coef = sprintf("spd%02d", 1), num = Inf)[
# tmp_genes,
# ]

gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
) |> select(
  ensembl, gene
)

gene_names <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    "Escher_genelist.xlsx"
  ),
  col_names = FALSE
) |> unlist()
names(gene_names) <- "gene"

# Each SPD ----

for (spd_it in colnames(cont.mat)) {
  # browser()
  spd_cont_df <- topTable(contrast_fit, coef = spd_it, num = Inf) |>
    mutate(sig_10 = adj.P.Val <= 0.1) |>
    rownames_to_column("gene_id") |>
    # Find gene name
    left_join(gene_df,
      by = c("gene_id" = "ensembl")
    )

  sig_gene_df <- spd_cont_df |> filter(gene %in% gene_names)


  # Make volcano plot ----
  pdf(
    here(
      "plots/PB_dx_genes/volcano_plot",
      paste0("volcano_plot_", spd_it, ".pdf")
    )
  )
  print(
    ggplot(
      data = spd_cont_df,
      aes(x = logFC, y = -log10(P.Value), color = sig_10)
    ) +
      geom_point() +
      geom_label_repel(
        data = sig_gene_df, # Add labels last to appear as the top layer
        aes(label = gene),
        force = 2,
        nudge_y = 0.1
      ) +
      # geom_label_repel(
      #   data = out_gene_df,
      #   aes(label = gene),
      #   color = "red",
      #   force = 2,
      #   nudge_y = -0.1
      # ) +
      labs(
        title = spd_it,
        y = "-log10(p-value)",
        x = "log2(FC)"
      ) +
      scale_color_manual(
        values = c(
          "FALSE" = "grey",
          "TRUE" = "red"
        ),
        breaks = c("TRUE", "FALSE"),
        name = "FDR < 0.1"
      ) +
      theme_minimal()
  )
  dev.off()
}
