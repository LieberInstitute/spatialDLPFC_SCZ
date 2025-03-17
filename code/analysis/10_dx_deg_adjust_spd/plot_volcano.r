suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(sessioninfo)
})

# Load Data ----
p_cut_off <- 0.10

gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd/preliminary",
    "test_PRECAST_07.csv"
  )
  # here(
  #   "processed-data/rds/10_dx_deg_adjust_spd",
  #   TODO: correct the file name
  #   "test_PRECAST_07.csv"
  # )
) |> mutate(
  sig_05 = fdr_scz <= 0.05,
  sig_10 = fdr_scz <= 0.10
)

gene_names <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    "Escher_genelist.xlsx"
  ),
  col_names = FALSE
) |> unlist()
names(gene_names) <- "gene"

sig_gene_df <- gene_df |> filter(gene %in% gene_names)


# sig_gene_df <- gene_df |>
#   filter(fdr_ntc <= p_cut_off)

# Make volcano plot ----
pdf(
  here(
    "plots/PB_dx_genes/",
    "volcano_plot_all.pdf"
  )
)
ggplot(
  data = gene_df,
  aes(x = logFC_scz, y = -log10(p_value_scz), color = sig_10)
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
    title = "Overall",
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
dev.off()


# Session Info ----
session_info()
