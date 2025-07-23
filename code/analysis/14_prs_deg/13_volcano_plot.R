# Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(sessioninfo)
  library(here)
})

# Load data ----
test_res <- read_csv(
  file = here(
    "processed-data/rds/14_prs_deg",
    "PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
) |>
  mutate(
    sig_10 = adj.P.Val < 0.10,
    gene_cat = case_when(
      sig_10 == FALSE ~ "grey",
      sig_10 == TRUE & logFC >= 0 ~ "red",
      sig_10 == TRUE & logFC < 0 ~ "blue"
    )
  )

## Create a volcano plot ----
ggplot(
  test_res,
  aes(x = logFC, y = -log10(P.Value))
) +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_point(
    aes(color = gene_cat)
  ) +
  scale_color_manual(
    values = c(
      "grey" = "grey",
      "red" = "red",
      "blue" = "blue"
    ),
    name = "Genes",
    labels = c(
      "grey" = "Not sig.",
      "red" = "Up (FDR < 0.10)",
      "blue" = "Down (FDR < 0.10)"
    )
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(
  filename = here(
    "plots/14_prs_deg",
    "PRS_DEG_volcano_plot_PRECAST07_donor_spd.pdf"
  ),
  plot = last_plot(),
  width = 4, height = 4.5
)

# Session Info ----
session_info()
