# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(ggrepel)
  library(sessioninfo)
})

# Load data ----
## Load SCZ-DEG ----
gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
) #|>
# Filter to FDR<0.10 genes
# filter(fdr_scz < 0.10)

# Error prevention
# stopifnot(nrow(gene_df) == 172)

## Load BrainSeq V2 ----
load(
  here(
    "processed-data/BrainSeqV2",
    "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda"
  )
)
# NOTE: contains outGene (what we need w. qSV) and outGene0 and outGeneNoAdjust
# https://github.com/LieberInstitute/brainseq_phase2?tab=readme-ov-file#dxstats_dlpfc_filtered_qsva_genelevel_nohgoldqsv_matchdlpfcrda

brainseq_v2_df <- outGene

nrow(brainseq_v2_df)
# [1] 24652
# NOTE: transcriptome wide

names(brainseq_v2_df)
#  [1] "Length"       "gencodeID"    "ensemblID"    "gene_type"
#  [5] "Symbol"       "EntrezID"     "Class"        "meanExprs"
#  [9] "NumTx"        "gencodeTx"    "passExprsCut" "logFC"
# [13] "AveExpr"      "t"            "P.Value"      "adj.P.Val"
# [17] "B"

hist(brainseq_v2_df$P.Value, breaks = 100)

sum(brainseq_v2_df$adj.P.Val < 0.1, na.rm = TRUE)
# [1] 632

### BrainSeq V2 statistics ----
tmp <- brainseq_v2_df |>
filter(adj.P.Val < 0.1)
# Number of down-reg sig genes.
sum(tmp$t < 0, na.rm = TRUE)


## Create merged data set ----
# Inner join data as experiment
inner_merge_df <- inner_join(
  gene_df |> filter(fdr_scz < 0.10),
  brainseq_v2_df |> filter(adj.P.Val < 0.10),
  by = c("ensembl" = "ensemblID")
)

nrow(inner_merge_df)
# [1] 36

# Full joined data
merged_df <- left_join(
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


hl_genes <- merged_df |>
  filter(Symbol %in% c("MAPK3", "ATP2B4", "KCNK1"))


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
      "Neither" = "lightgrey"
    )
  ) +
  labs(
    x = "SCZ-DEG T-statistic",
    y = "BrainSeq V2 T-statistic",
    title = sprintf("Scatter plot of T-statistics (N= %d genes)", nrow(merged_df)),
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
  filter(adj.P.Val < 0.1 & fdr_scz < 0.10) |>
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

# Venn Diagram ----
library(VennDiagram)

# Create a Venn diagram
venn.diagram(
  x = list(
    `172 DEGs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `BrainSeq V2` = brainseq_v2_df |> pull(ensemblID)
  ),
  filename = here("plots/10_dx_deg_adjust_spd/venn_diagram_compare_brainseqv2.tiff"),
  disable.logging = TRUE,
  fill = c("orange", "purple"),
  alpha = 0.5,
  cex = 1,
  cat.cex = 1,
  cat.pos = 0,
  height = 1.5,
  width = 3,
  units = "in",
  resolution = 300
)

# SynGO analysis ----
## Load SynGO results ----
synGO_df <- read_xlsx(
  here(
    "processed-data/SynGO/syngo_genes.xlsx"
  )
)

merged_df <- merged_df |>
  mutate(
    synGO = ensembl %in% synGO_df$ensembl_id
  )

# Prob that 172 DEG is a synGO
merged_df |>
  filter(fdr_scz < 0.10) |>
  pull(synGO) |> 
  mean()
# [1] 0.127907

# Prob that BrainSeqV2 is a synGO
merged_df |>
  filter(adj.P.Val < 0.10) |>
  pull(synGO) |> 
  mean()
# [1] 0.1585624

# Odds ratio that 
(0.1586/(1-0.1586))/(0.1279/(1-0.1279))
# [1] 1.285276


merged_df |>
  filter(adj.P.Val > 0.1 & fdr_scz < 0.10) |> 
  write_csv(
    here(
      "code/analysis/10_dx_deg_adjust_spd",
      "dx-DEGs_on_sig_in_PNN.csv"
    )
  )

merged_df |>
  filter(adj.P.Val < 0.1 & fdr_scz > 0.10) |> 
  write_csv(
    here(
      "code/analysis/10_dx_deg_adjust_spd",
      "dx-DEGs_on_sig_in_BrainSeqV2.csv"
    )
  )




# Session Info ----
sessioninfo::session_info()
