# Load library  ----
suppressForeignCheck({
  library(here)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

# Load Data ----
## Load DEG results ----
dx_res <- read_csv(
  here(
    "code/analysis/dx_deg_spg_neuropil",
    "neuropil-dx_DEG-GM.csv"
  )
) |>
  mutate(
    sig_05 = fdr_scz <= 0.05,
    sig_10 = fdr_scz <= 0.10,
    gene_cat = case_when(
      sig_05 == FALSE ~ "grey",
      sig_05 == TRUE & logFC_scz >= 0 ~ "red",
      sig_05 == TRUE & logFC_scz < 0 ~ "blue"
    )
  ) |>
  arrange(
    desc(p_value_scz)
  )

## Label data ----
sig_gene_df <- dx_res |>
  filter(gene %in%
    c(
      "MAPK3", "PFN1", "RTN4", "VGF",
      "CORT", "SST"
    ))

# Structure Volcano plots ----
ggplot(
  data = dx_res,
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
    data = sig_gene_df, # Add labels last to appear as the top layer
    aes(x = logFC_scz, y = -log10(p_value_scz), label = gene, color = gene_cat),
    force = 2,
    nudge_y = 0.1,
    size = 3
  ) +
  # Format the legends
  scale_color_identity(
    name = "Significance",
    labels = c(
      "FDR > 0.10",
      "Up-regulated",
      "Down-regulated"
    )
  ) +
  labs(
    title = sprintf(
      "SCZ-Differentially Expressed Genes among Neuropil+ Spots (N = %d)",
      nrow(dx_res)
    ),
    y = "-log10(p-value)",
    x = "log2(Fold Change)",
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

## Save Plot ----
ggsave(
  here(
    "code/analysis/dx_deg_spg_neuropil",
    "volcano_plot_neuropil-dx_DEG-GM.pdf"
  )
)


# Session Info ----
session_info()
