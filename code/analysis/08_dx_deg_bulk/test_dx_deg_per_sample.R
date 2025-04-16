# Load Libray --------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(sessioninfo)
  library(tidyverse)
  library(ggrepel)
})

# Load data -----

spe_bulk <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_donor.rds"
  )
)

# TODO: find the pseudobulk data

# Pseudo-bulk DE ----
spe_bulk <-
  registration_pseudobulk(
    spe,
    var_registration = "dx",
    var_sample_id = "sample_id",
    min_ncells = 10
  )

covars <- c("age", "sex", "lot_num")

# set.seed(20240411)
# spe_bulk <- runPCA(sce_pseudo)

# tmp <- registration_stats_enrichment(
#   spe_bulk,
#   block_cor = NaN,
#   covars = c("age", "sex"),
#   # var_registration = "dx",
#   gene_ensembl = "gene_id",
#   gene_name = "gene_name"
# )

# tmp |>
#   filter(fdr_scz <= 0.05) |>
#   nrow()

# tmp <- registration_stats_enrichment(
#   spe_bulk,
#   block_cor = NaN,
#   covars = c("age", "sex", "lot_num"),
#   # var_registration = "dx",
#   gene_ensembl = "gene_id",
#   gene_name = "gene_name"
# )



# plotPCA(
#   spe_bulk,
#   colour_by = "spd",
#   ncomponents = 6,
#   point_size = 0.5,
#   label_format = c("%s %02i", " (%i%%)"),
# )





registration_mod <-
  registration_model(spe_bulk,
    covars = covars
  )


block_cor <-
  registration_block_cor(
    spe_bulk,
    registration_model = registration_mod
  )

results_enrichment <-
  registration_stats_enrichment(
    spe_bulk,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )

gene_df <- results_enrichment

nrow(gene_df)

saveRDS(
  results_enrichment,
  file = here(
    # here
    # "processed-data/rds/pseudo_bulk/sample_bulk_DE.rds"
    )
)

# Bench mark to BrainSeqV2
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


merged_df <- inner_join(
  gene_df,
  brainseq_v2_df,
  by = c("ensembl" = "ensemblID")
) |>
  mutate(
    study_sig = case_when(
      fdr_SCZ > 0.1 & adj.P.Val < 0.1 ~ "BrainSeq V2",
      fdr_SCZ < 0.1 & adj.P.Val > 0.1 ~ "SCZ-DEG",
      fdr_SCZ < 0.1 & adj.P.Val < 0.1 ~ "Both",
      TRUE ~ "Neither"
    )
  )


hl_genes <- merged_df |>
  filter(Symbol %in% c("MAPK3", "ATP2B4", "KCNK1"))


# Make scatter plots ----
## Scatter plot of sig genes ----
ggplot(merged_df, aes(x = t_stat_SCZ, y = t, color = study_sig)) +
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
        merged_df$t_stat_SCZ,
        merged_df$t,
        use = "pairwise.complete.obs"
      )
    ),
    hjust = -0.1, vjust = 2, size = 5, color = "black"
  )

# Comparing to analysis accounting for spatial domain ----
# Load results ----
dx_res_spd <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)

merged_df_spd <- inner_join(
  gene_df,
  dx_res_spd,
  by = c("ensembl"),
  suffix = c("_donor", "_spd")
) #|>
  # mutate(
  #   study_sig = case_when(
  #     fdr_SCZ > 0.1 & adj.P.Val < 0.1 ~ "BrainSeq V2",
  #     fdr_SCZ < 0.1 & adj.P.Val > 0.1 ~ "SCZ-DEG",
  #     fdr_SCZ < 0.1 & adj.P.Val < 0.1 ~ "Both",
  #     TRUE ~ "Neither"
  #   )
  # )

ggplot(merged_df_spd, aes(x = t_stat_SCZ, y = t_stat_scz)) +
  # Prepare grid
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
  # Add genes
  geom_point(alpha = 0.7) +
  # geom_label_repel(
  #   data = hl_genes,
  #   aes(label = Symbol),
  #   box.padding = 0.5,
  #   point.padding = 0.5,
  #   segment.color = "grey50",
  #   segment.size = 0.5,
  #   segment.alpha = 0.5,
  #   size = 3
  # ) +
  # # Format plot
  # scale_color_manual(
  #   values = c(
  #     "Both" = "#D55E00",
  #     "BrainSeq V2" = "#0072B2",
  #     "SCZ-DEG" = "#009E73",
  #     "Neither" = "lightgrey"
  #   )
  # ) +
  labs(
    x = "donor_spd T-statistic",
    y = "Donor_level T-statistic",
    title = sprintf("Scatter plot of T-statistics (N= %d genes)", nrow(merged_df_spd)),
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
        merged_df_spd$t_stat_SCZ,
        merged_df_spd$t_stat_scz,
        use = "pairwise.complete.obs"
      )
    ),
    hjust = -0.1, vjust = 2, size = 5, color = "black"
  )

# Comparing to analysis accounting for spatial domain ----
# Load results ----
dx_res_spd <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)



# Session Info ----
sessioninfo::session_info()
