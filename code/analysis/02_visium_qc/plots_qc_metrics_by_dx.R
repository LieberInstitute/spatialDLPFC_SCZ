# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(scater)
  library(ggpubr)
  library(sessioninfo)
})

# Load Data ----
## Load Spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "raw_spe_wo_SPG_N63.rds"
  )
)

## Add sample info ----
spe$DX <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]
spe$DX <- toupper(spe$DX)

spe$brnum <- metadata(spe)$dx_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$DX)
)

stopifnot(
  length(unique(spe$sample_label)) == 63
)

## Load outlier keys ----
tot_outlier_df <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "combined_outlier_df.rds"
  )
)

## Merge data together ----
spe$all_outlier <- tot_outlier_df[spe$key, "all_outlier"]
spe$discard <- tot_outlier_df[spe$key, "remove"]

## Remove out_tissue_spots ----
spe <- spe[, spe$in_tissue == TRUE]

# Plots ----
# NOTES: Becasue of the huge amount of datapoints,
# save plots to png saves a lot of time to view them
# during the development and as a finalized figure.

## sum_umi ----
p_umi <- plotColData(spe,
  x = "sample_label", y = "sum_umi",
  color = "discard",
  point_size = 0.1,
  scattermore = TRUE
) +
  geom_hline(
    aes(yintercept = 100),
    alpha = 0.8,
    linetype = "dashed", color = "lightgrey"
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~ spe$DX, scales = "free_x") +
  labs(
    title = "Total UMI",
    y = "Library Size",
    x = "Sample"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.justification = "center",
    legend.position = "none"
  ) +
  scale_color_manual(
    name = "Spots",
    limits = c("TRUE", "FALSE"), # Show discard first
    values = c(
      "TRUE" = "deeppink",
      "FALSE" = "grey90"
    ),
    labels = c(
      "TRUE" = "Discarded",
      "FALSE" = "Kept"
    )
  ) +
  # Make legend circle more obvious
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(
  "test_plot.png",
  p_umi,
  width = 9, height = 4, units = "in"
)

## sum_gene ----
p_gene <- plotColData(spe,
  x = "sample_label", y = "sum_gene",
  color = "discard",
  point_size = 0.1,
  scattermore = TRUE
) +
  geom_hline(
    aes(yintercept = 200),
    alpha = 0.8,
    linetype = "dashed", color = "lightgrey"
  ) +
  facet_wrap(~ spe$DX, scales = "free_x") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    title = "Unique Genes",
    y = "Number of Genes",
    x = "Sample"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.justification = "center",
    legend.position = "none"
  ) +
  scale_color_manual(
    name = "Spots",
    limits = c("TRUE", "FALSE"), # Show discard first
    values = c(
      "TRUE" = "deeppink",
      "FALSE" = "grey90"
    ),
    labels = c(
      "TRUE" = "Discarded",
      "FALSE" = "Kept"
    )
  ) +
  # Make legend circle more obvious
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(
  "test_plot.png",
  p_gene,
  width = 9, height = 4, units = "in"
)

## Mito_ratio ----
p_mito <- plotColData(spe,
  x = "sample_label", y = "expr_chrM_ratio",
  color = "discard",
  point_size = 0.1,
  scattermore = TRUE
) +
  facet_wrap(~ spe$DX, scales = "free_x") +
  labs(
    title = "Mito. Ratio",
    y = "Proportion",
    x = "Sample"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = 12, angle = 90, hjust = 1
    ),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.justification = "center",
    legend.position = "none"
  ) +
  scale_color_manual(
    name = "Spots",
    limits = c("TRUE", "FALSE"), # Show discard first
    values = c(
      "TRUE" = "deeppink",
      "FALSE" = "grey90"
    ),
    labels = c(
      "TRUE" = "Discarded",
      "FALSE" = "Kept"
    )
  ) +
  # Make legend circle more obvious
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(
  "test_plot.png",
  p_mito,
  width = 9, height = 4, units = "in"
)

## Create panel ----
paneled_p <- ggpubr::ggarrange(
  p_umi, p_gene, p_mito,
  ncol = 1,
  heights = c(2, 2, 2.5),
  align = "v",
  legend = "bottom",
  common.legend = TRUE
)

## Save plot ----
ggsave(
  here(
    "plots/02_visium_qc",
    "qc_plot_metrics_distribution.png"
  ),
  plot = paneled_p,
  width = 8, height = 10, unit = "in"
)

# Session Info ----
sessioninfo::session_info()
