# Load library -----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggrepel)
  library(ggpubr)
  library(sessioninfo)
})

# Load DEG result ----
## Load layer-adjusted dx-DEGs ----
dx_deg_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)

## Load layer-adjusted PRS-DEGs ----
prs_deg_df <- read_csv(
  file = here(
    # "processed-data/rds/14_prs_deg",
    # "PRS_DEG_test_res_PRECAST07_donor_spd.csv"
    "processed-data/rds/14_prs_deg",
    "norm_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)


## Merge results ----
merged_df <- dx_deg_df |>
  inner_join(
    prs_deg_df |>
      select(
        ensembl = gene_id, logFC_prs = logFC, AveExpr_prs = AveExpr,
        t_prs = t, P.Value_prs = P.Value, fdr_prs = adj.P.Val
      ),
    by = "ensembl",
    relationship = "one-to-one"
  ) |>
  mutate(
    study_sig = case_when(
      fdr_scz > 0.1 & fdr_prs < 0.1 ~ "PRS-DEG only",
      fdr_scz < 0.1 & fdr_prs > 0.1 ~ "DX-DEG only",
      fdr_scz < 0.1 & fdr_prs < 0.1 ~ "Both",
      TRUE ~ "Neither"
    ) |> factor(levels = c("Neither", "DX-DEG only", "PRS-DEG only", "Both"))
  )


### Descriptive statistics ----
cor(merged_df$t_stat_scz, merged_df$t_prs)
# [1] 0.5551971


# Scatter plot ----
hl_genes <- c(
  # "EIF2A", "COX7A2", "MCL1", "AIF1", "TREM2",
  # "VEGFA", "A2M", "FOSB", "JUN", "SYT1", "KANSL1-AS1",
  # "MAP2", "GAP43", "MAPK3", "MSI2", "R3HDM2", "THOC7",
  # "C1QA", "C1QB", "C3", "CX3CR1", "TYROBP", "CD74"
  "MAPK3",
  "AIF1", "TREM2", "C1QA", "C1QB",
   "C3", "CX3CR1", "TYROBP", "CD74", "KANSL1-AS1"
)

ret_p <- ggplot(
  merged_df |> arrange(study_sig),
  aes(x = t_stat_scz, y = t_prs)
) +
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(aes(color = study_sig)) +
  geom_text_repel(
    data = merged_df |>
      filter(gene %in% hl_genes),
    aes(label = gene),
    color = "black",
    force = 4,
    # All labels not overlapping
    max.overlaps = Inf,
    # Have arrows
    min.segment.length = 0,
    size = 2, # 6pt label font ≈ size 2 in ggplot2
    arrow = arrow(length = unit(2, "points"), type = "closed"),
    segment.size = 0.2,
    segment.color = "black"
  ) +
  annotate("text",
    x = 4, y = -5,
    label = sprintf(
      "r=%.3f",
      # Correlation of the t-statistics
      cor(
        merged_df$t_stat_scz, merged_df$t_prs,
        use = "pairwise.complete.obs"
      )
    ),
    hjust = -0.1, vjust = 2, size = 3, color = "black"
  ) +
  labs(
    x = "t-statistics (DX-DEGs)",
    y = "t-statistics (PRS-DEGs)"
  ) +
  scale_x_continuous(limits = c(-7, 7), breaks = seq(-6, 6, by = 2)) +
  scale_y_continuous(limits = c(-7, 7), breaks = seq(-6, 6, by = 2)) +
  scale_color_manual(
    values = c(
      "Neither" = "lightgrey",
      "DX-DEG only" = "#009E73",
      "PRS-DEG only" = "#0072B2",
      "Both" = "#D55E00"
    )
  ) +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.text = element_text(size = 6)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, label.position = "right", label.hjust = 0))



ggsave(
  filename = here(
    "plots/14_prs_deg",
    "cor_plot_norm_PRS_vs_Dx-DEGs_PRECAST07_donor_spd.pdf"
  ),
  plot = ret_p,
  width = 2.5, height = 2.7, units = "in"
)

# Session info ----
session_info()
