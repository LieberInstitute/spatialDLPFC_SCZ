# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(here)
  library(scater)
  library(escheR)
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

## Remove out_tissue_spots ---
spe <- spe[, spe$in_tissue == TRUE]

## Calcualte normaled counts ----
spe <- spe[, colSums(counts(spe)) != 0]
spe <- logNormCounts(spe)

## Adjust sample name to Brnum_dx ---
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$brnum <- metadata(spe)$dx_df$subject[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sample_label <- paste0(
  spe$brnum, "_", toupper(spe$dx)
)

## Subset to representative samples ----
# error prevention
stopifnot(all(c("V13M06-342_D1", "V13M06-343_D1") %in% unique(spe$sample_id)))

spe <- spe[, spe$sample_id %in% c("V13M06-342_D1", "V13M06-343_D1")]

# Create visualization ----
plot_list <- unique(spe$sample_label) |>
  set_names() |>
  map(.f = function(.smp) {
    sub_spe <- spe[, spe$sample_label == .smp]

    RBG_norm_MBP <- scales::rescale(
      logcounts(sub_spe)[which(rowData(sub_spe)$gene_name == "MBP"), ] |>
        scale(center = TRUE, scale = FALSE) |>
        sapply(max, 0),
      to = c(0, 1)
    )

    RBG_norm_PCP4 <- scales::rescale(
      logcounts(sub_spe)[which(rowData(sub_spe)$gene_name == "PCP4"), ] |>
        scale(center = TRUE, scale = FALSE) |>
        sapply(max, 0),
      to = c(0, 1)
    )

    RBG_norm_SNAP25 <- scales::rescale(
      logcounts(sub_spe)[which(rowData(sub_spe)$gene_name == "SNAP25"), ] |>
        scale(center = TRUE, scale = FALSE) |>
        sapply(max, 0), # ignore small counts
      to = c(0, 1)
    )

    sub_spe$rbg_val <- rgb(
      red = RBG_norm_MBP,
      green = RBG_norm_PCP4,
      blue = RBG_norm_SNAP25
    )

    # browser()
    make_escheR(sub_spe) |>
      add_fill(
        "rbg_val",
        point_size = 2.1
      ) +
      scale_fill_identity() +
      labs(title = .smp) +
      theme(
        # legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )
  })


# Create a custom legend ----
legend_plot <- data.frame(
  x = 1, y = 1,
  Gene = c("MBP", "PCP4", "SNAP25")
) |>
  ggplot() +
  geom_point(aes(x, y, color = Gene)) + # Increase the size of the points
  scale_color_manual(
    name = "Gene",
    values = c(
      "MBP" = "red",
      "PCP4" = "green",
      "SNAP25" = "blue"
    ),
    guide = guide_legend(
      override.aes = list(size = 7) # Increase the size of the legend keys
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )

# Combine plots with legend ----
combined_plot <- ggpubr::ggarrange(
  plotlist = plot_list,
  nrow = 1,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = ggpubr::get_legend(legend_plot)
)

# Save combined plot ----
ggsave(
  filename = here(
    "plots/02_visium_qc",
    "qc_spot_plot_anatomical_orientation_rep_sample_with_legend.pdf"
  ),
  plot = combined_plot,
  height = 5.5, width = 9.5,
  unit = "in"
)

# Session Info
session_info()
