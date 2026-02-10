# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(sessioninfo)
  library(ggrepel)
  library(here)
})

# Load data ----
test_res <- read_csv(
  file = here(
    "processed-data/rds/14_prs_deg",
    # "PRS_DEG_test_res_PRECAST07_donor_spd.csv"
    # "processed-data/rds/14_prs_deg/norm_PRS_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
) |>
  mutate(
    sig_10 = adj.P.Val < 0.10,
    gene_cat = case_when(
      sig_10 == FALSE & P.Value < 0.05 & logFC >= 0 ~ "#ef8080",
      sig_10 == FALSE & P.Value < 0.05 & logFC < 0 ~ "#808ef8",
      sig_10 == FALSE ~ "grey",
      sig_10 == TRUE & logFC >= 0 ~ "red",
      sig_10 == TRUE & logFC < 0 ~ "blue"
    )
  )

## Create a volcano plot ----
hl_genes <- c(
  "EIF2A", "COX7A2", "MCL1", "AIF1", "TREM2",
  "VEGFA", "A2M", "FOSB", "JUN", "SYT1", "KANSL1-AS1",
  "MAP2", "GAP43", "MAPK3", "MSI2", "R3HDM2", "THOC7"
)
length(hl_genes)

ggplot(
  test_res |> arrange(desc(P.Value)),
  aes(x = logFC, y = -log10(P.Value))
) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey"
  ) +
  geom_point(
    aes(color = gene_cat)
  ) +
  geom_text_repel(
    data = test_res |>
      filter(gene_name %in% hl_genes),
    aes(label = gene_name),
    color = "black",
    force = 0.5,
    # All labels not overlapping
    max.overlaps = Inf,
    # Have arrows
    min.segment.length = 0,
    size = 2, # 6pt label font ≈ size 2 in ggplot2
    arrow = arrow(length = unit(0.5, "points"), type = "closed"),
    segment.size = 0.2,
    segment.color = "black"
  ) +
  # scale_color_manual(
  #   values = c(
  #     "grey" = "grey",
  #     "red" = "red",
  #     "blue" = "blue"
  #   ),
  #   name = "Genes",
  #   labels = c(
  #     "grey" = "Not sig.",
  #     "red" = "Up (FDR < 0.10)",
  #     "blue" = "Down (FDR < 0.10)"
  #   )
  # ) +
  theme_classic(base_size = 6) +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "bottom"
  ) +
  scale_color_identity(
    # name = "Significance",
    # labels = c(
    #   "FDR > 0.10",
    #   "Up-regulated",
    #   "Down-regulated"
    # ),
    # guide = "legend"
  ) +
  labs(
    # title = "PRS-associated DEGs",
    # y = "-log10(p-value)",
    x = "log2(FC in PRS)"
  )

ggsave(
  filename = here(
    # "plots/14_prs_deg",
    # "PRS_DEG_volcano_plot_PRECAST07_donor_spd.pdf"
    "plots/14_prs_deg",
    "norm_PRS_DEG_volcano_plot_PRECAST07_donor_spd.pdf"
  ),
  plot = last_plot(),
  width = 2.5, height = 2.5, units = "in"
)

# Session Info ----
session_info()
