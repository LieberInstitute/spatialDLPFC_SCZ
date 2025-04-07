# Load library  ----
suppressForeignCheck({
  library(here)
  library(spatialLIBD)
  library(limma)
  library(readxl)
  library(ggrepel)
  library(sessioninfo)
})


# Load data ----
dx_res <- read_csv(
  here(
    "code/analysis/dx_deg_spg_neuropil",
    "neuropil-dx_DEG-GM.csv"
  )
) |>
  mutate(
    gene_cat = case_when(
      fdr_scz <= 0.05 & logFC_scz > 0 ~ "red",
      fdr_scz <= 0.05 & logFC_scz < 0 ~ "blue",
      TRUE ~ "grey"
    )
  )


synGO_genes <- read_excel(
  here::here(
    "processed-data/SynGO/syngo_genes.xlsx"
  )
)


sig_synGO_gene_df <- dx_res |>
  filter(
    ensembl %in% synGO_genes$ensembl_id & fdr_scz <= 0.05
  )



ggplot(
  data = dx_res |>
    mutate(),
  aes(x = logFC_scz, y = -log10(p_value_scz))
) +
  # Add line for nominal threshold 0.05
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "darkgrey"
  ) +
  # Plot all genes
  geom_point(
    aes(color = gene_cat),
    size = 2
  ) +
  # Add interested gene labels
  geom_label_repel(
    data = sig_synGO_gene_df, # Add labels last to appear as the top layer
    aes(
      x = logFC_scz, y = -log10(p_value_scz), label = gene, color = gene_cat
    ),
    force = 2,
    nudge_y = 0.1,
    size = 3
  ) +
  # Format the legends
  scale_color_identity(
    name = "Significance",
    labels = c(
      "Down-regulated",
      "FDR > 0.05",
      "Up-regulated"
    ),
    guide = "legend"
  ) +
  labs(
    title = "SCZ-Differentially Expressed Genes among Neuropil+",
    y = "-log10(p-value)",
    x = "log2(Fold Change in SCZ)",
    color = "Significance"
  ) +
  # Format the plot
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    panel.border = element_rect(
      color = "black", fill = NA, size = 1
    )
  )

ggsave(
  filename = here(
  "code/analysis/dx_deg_spg_neuropil",
  "neuropil_dx_DEG_GM.pdf"
  ),
  width = 8,
  height = 6
)
