# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

# Load data ----
## Load SCZ-DEG ----
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
) #|>
# Filter to FDR<0.10 genes
# filter(fdr_scz < 0.10)

# Error prevention
stopifnot(nrow(gene_df) == 172)

## Load BrainSeq V2 ----
brainseq_v2_df <- read_csv(
  here(
    "code/analysis/10_dx_deg_adjust_spd/",
    "Collado-Torre_2019-SupplementaryTable11_sczd_gene.csv"
  )
) |> filter(
  region == "DLPFC"
)

summary(brainseq_v2_df$P.Value)

range(brainseq_v2_df$adj.P.Val)


nrow(brainseq_v2_df)
# [1] 632

## Create merged data set ----
# Inner join data as experiment
inner_merge_df <- inner_join(
  gene_df,
  brainseq_v2_df,
  by = c("ensembl" = "ensemblID")
)

nrow(inner_merge_df)
# [1] 42

inner_merge_df |>
  select(
    ensembl, gene, gene_type, fdr_scz, adj.P.Val
  ) |>
  View()


# Full joined data
merged_df <- full_join(
  gene_df,
  brainseq_v2_df,
  by = c("ensembl" = "ensemblID")
) |>
  mutate(
    study_sig = case_when(
      fdr_scz > 0.1 & adj.P.Val < 0.1 ~ "BrainSeq V2",
      fdr_scz < 0.1 & adj.P.Val > 0.1 ~ "SCZ-DEG",
      fdr_scz < 0.1 & adj.P.Val < 0.1 ~ "Both",
      TRUE ~ "Neither"
    )
  )
stopifnot(!"HIPPO" %in% unique(merged_df$region))


hl_genes <- merged_df |>
  filter(Symbol %in% c("MAPK3"))


# Make scatter plots ----
## Scatter plot of sig genes ----
ggplot(merged_df, aes(x = t_stat_scz, y = t, color = study_sig)) +
  # Prepare grid
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
  # Add genes
  geom_point(alpha = 0.7) +
  geom_label_repel(
    data = hl_genes,
    aes(label = Symbol),
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey50",
    segment.size = 0.5,
    segment.alpha = 0.5,
    size = 3
  ) +
  # Format plot
  scale_color_manual(
    values = c(
      "Both" = "#D55E00",
      "BrainSeq V2" = "#0072B2",
      "SCZ-DEG" = "#009E73",
      "Neither" = "#F0E442"
    )
  ) +
  labs(
    x = "SCZ-DEG T-statistic",
    y = "BrainSeq V2 T-statistic",
    title = "Scatter plot of T-statistics",
    color = "Study Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  annotate("text",
    x = -Inf, y = Inf,
    label = sprintf(
      "r=%.3f",
      # Correlation of the t-statistics
      cor(
        merged_df$t_stat_scz,
        merged_df$t,
        use = "pairwise.complete.obs"
      )
    ),
    hjust = -0.1, vjust = 2, size = 5, color = "black"
  )

ggsave(
  here(
    "plots/10_dx_deg_adjust_spd",
    "scatter_plot_t_stat_brain_seq_v2.pdf"
  ),
  width = 8, height = 6
)


# Save plots ----


# Gene significant in both studies
merged_df |> dim()

# merged_df |> filter(adj.P.Val < 0.10)

merged_df |>
  filter(adj.P.Val & fdr_scz < 0.10) |>
  write.csv(
    here(
      "processed-data/rds/10_dx_deg_adjust_spd",
      "brainseq-v2_overlap_genes.csv"
    )
  )

# What are the genes sig in our study
merged_df |>
  filter(fdr_scz < 0.10 & is.na(adj.P.Val)) |>
  pull(gene)

# Session Info ----
sessioninfo::session_info()
