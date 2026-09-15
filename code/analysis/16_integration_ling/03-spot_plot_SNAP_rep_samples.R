# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(escheR)
  library(Polychrome)
  library(ggpubr)
  library(sessioninfo)
})


# Config ----
fld_plot <- here("plots/16_integration_ling")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)

# Representative samples: Br8667 (NTC) and Br5973 (SCZ)
# (same pair used throughout the paper, e.g. 03_visium_spatial_clustering)
rep_sample_id <- c("V13M06-342_D1", "V13M06-343_D1")


# Load Data ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

score_df <- readRDS(
  here(
    "processed-data/rds/16_integration_ling",
    "SNAP_score_df.rds"
  )
)

## Merge normalized SNAP scores into colData(spe) ----
col_data_df <- score_df |>
  select(key, SNAPa_zscore, SNAPn_zscore) |>
  right_join(
    colData(spe) |> data.frame(),
    by = "key",
    relationship = "one-to-one"
  )
rownames(col_data_df) <- col_data_df$key

# error prevention
stopifnot(identical(col_data_df$key, spe$key))

spe$SNAPa_zscore <- col_data_df$SNAPa_zscore
spe$SNAPn_zscore <- col_data_df$SNAPn_zscore

## sample_label used for panel titles, e.g. Br8667_NTC ----
spe$sample_label <- paste0(spe$brnum, "_", toupper(spe$dx))


# Subset to representative samples ----
# error prevention
stopifnot(all(rep_sample_id %in% unique(spe$sample_id)))

spe <- spe[, spe$sample_id %in% rep_sample_id]

# error prevention
stopifnot(length(unique(spe$sample_id)) == 2)

# error prevention
stopifnot("factor" %in% class(spe$spd_label))


# Spatial domain color palette ----
spd_palette <- set_names(
  Polychrome::palette36.colors(length(levels(spe$spd_label))),
  levels(spe$spd_label)
)


# Make spot plot ----
## Function to create escheR spot plot for a given normalized SNAP score ----
plot_snap_spot <- function(spe, score_var, title) {
  plot_list <- unique(spe$sample_label) |>
    set_names() |>
    map(.f = function(.smp) {
      sub_spe <- spe[, spe$sample_label == .smp]

      make_escheR(sub_spe) |>
        add_fill(
          score_var,
          point_size = 2.1
        ) |>
        add_ground(
          "spd_label",
          point_size = 2.2
        ) +
        labs(title = .smp) +
        theme(
          plot.title = element_text(size = 20, hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)
        )
    })

  ## Shared, diverging color scale (score is z-scored, centered at 0) ----
  score_range <- range(spe[[score_var]], na.rm = TRUE)
  score_lim <- max(abs(score_range)) * c(-1, 1)

  plot_list <- plot_list |> lapply(FUN = function(.p) {
    .p +
      scale_fill_gradient2(
        name = title,
        limits = score_lim,
        low = "steelblue", mid = "white", high = "firebrick",
        midpoint = 0
      ) +
      scale_color_manual(
        name = "Spatial Domain",
        values = spd_palette,
        guide = guide_legend(override.aes = list(size = 7))
      ) +
      theme(
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
      )
  })

  ## Combine plots with legend ----
  # Both panels share the same fill limits and spatial domain palette, so a
  # single, common legend applies to both.
  ggpubr::ggarrange(
    plotlist = plot_list,
    nrow = 1,
    ncol = 2,
    common.legend = TRUE,
    legend = "right"
  )
}


## SNAP-a spot plot ----
combined_plot_snapa <- plot_snap_spot(
  spe,
  score_var = "SNAPa_zscore",
  title = "SNAP-a\n(normalized)"
)

ggsave(
  filename = file.path(fld_plot, "spot_plot_SNAPa_rep_samples.pdf"),
  plot = combined_plot_snapa,
  height = 5.5, width = 11.5,
  units = "in"
)

## SNAP-n spot plot ----
combined_plot_snapn <- plot_snap_spot(
  spe,
  score_var = "SNAPn_zscore",
  title = "SNAP-n\n(normalized)"
)

ggsave(
  filename = file.path(fld_plot, "spot_plot_SNAPn_rep_samples.pdf"),
  plot = combined_plot_snapn,
  height = 5.5, width = 11.5,
  units = "in"
)

print("Finished making SNAP-a / SNAP-n spot plots for representative samples")


# Session Info ----
sessioninfo::session_info()
