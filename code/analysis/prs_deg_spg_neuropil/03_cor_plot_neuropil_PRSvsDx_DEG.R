# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(sessioninfo)
})

# Load Data ----
## Neuropil Dx-DEG results ----
dx_deg_df <- read_csv(
  here(
    "code/analysis/dx_deg_spg_neun",
    "neun-dx_DEG-GM.csv"
  )
)

## Neuropil PRS-DEG result ----
prs_deg_df <- read_csv(
  here(
    "processed-data/rds/14_prs_deg",
    "Neuropil_PRS_DEG_test_res_PRECAST07_donor_spd.csv"
  )
)

## Merged df -----
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
# [1] 0.3806872


# Scatter plot ----
ggplot(
  merged_df |> arrange(study_sig),
  aes(x = t_stat_scz, y = t_prs)
) +
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(aes(color = study_sig)) +
  geom_text_repel(
    data = merged_df |>
      filter(gene %in% c(
        "XRRA1", "C2orf74", "AC004556.3", "MTRNR2L8", # PRS-DEG (all),
        "SURF1", "MRPL23" # PRS-DEG in Neuropil
        # "MAPK3", "C3", "AIF1", "SERPINA3", "A2M", # BrainseV2 & dx-DEG
        # "BDNF", "FKBP5", "SST", "GFAP", "MT2A" # dx-DEG
      )),
    aes(label = gene),
    color = "black",
    force = 4,
    # All labels not overlapping
    max.overlaps = Inf,
    # Have arrows
    min.segment.length = 0,
    size = 4, # 6pt label font ≈ size 2 in ggplot2
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
    hjust = -0.1, vjust = 2, size = 4, color = "black"
  ) +
  labs(
    title = "Correlation of t-statistics: Neuropil PRS-DEGs vs Dx-DEGs",
    x = "t-statistics (DX-DEGs among Neuropil+)",
    y = "t-statistics (PRS-DEGs among Neuropil+)"
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
  theme_classic(base_size = 12)

ggsave(
  filename = here(
    "plots/14_prs_deg",
    "cor_plot_Neuropil_PRS_vs_Dx-DEGs.pdf"
  ),
  width = 7, height = 5.5, units = "in"
)

# Session info ----
session_info()
