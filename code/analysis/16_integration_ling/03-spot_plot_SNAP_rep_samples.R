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

## Merge SNAP scores into colData(spe) ----
col_data_df <- score_df |>
  select(
    key,
    SNAPa_top_2000_score, SNAPn_top_2000_score,
    SNAPa_all_genes_score, SNAPn_all_genes_score
  ) |>
  right_join(
    colData(spe) |> data.frame(),
    by = "key",
    relationship = "one-to-one"
  )
rownames(col_data_df) <- col_data_df$key

# error prevention
stopifnot(identical(col_data_df$key, spe$key))

spe$SNAPa_top_2000_score <- col_data_df$SNAPa_top_2000_score
spe$SNAPn_top_2000_score <- col_data_df$SNAPn_top_2000_score
spe$SNAPa_all_genes_score <- col_data_df$SNAPa_all_genes_score
spe$SNAPn_all_genes_score <- col_data_df$SNAPn_all_genes_score

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
## Function to create escheR spot plot for a given SNAP score ----
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

  ## Shared color scale (grayscale: white = min, black = max) ----
  score_range <- range(spe[[score_var]], na.rm = TRUE)

  fill_scale <- scale_fill_gradient(
    name = title,
    limits = score_range,
    low = "white", high = "black"
  )

  plot_list <- plot_list |> lapply(FUN = function(.p) {
    .p +
      fill_scale +
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


## SNAP-a top-2000 spot plot ----
combined_plot_snapa_top2000 <- plot_snap_spot(
  spe,
  score_var = "SNAPa_top_2000_score",
  title = "SNAP-a\n(top 2000)"
)

ggsave(
  filename = file.path(fld_plot, "spot_plot_SNAPa_top_2000_rep_samples.pdf"),
  plot = combined_plot_snapa_top2000,
  height = 5.5, width = 11.5,
  units = "in"
)

## SNAP-n top-2000 spot plot ----
combined_plot_snapn_top2000 <- plot_snap_spot(
  spe,
  score_var = "SNAPn_top_2000_score",
  title = "SNAP-n\n(top 2000)"
)

ggsave(
  filename = file.path(fld_plot, "spot_plot_SNAPn_top_2000_rep_samples.pdf"),
  plot = combined_plot_snapn_top2000,
  height = 5.5, width = 11.5,
  units = "in"
)

## SNAP-a all-genes spot plot ----
combined_plot_snapa_all_genes <- plot_snap_spot(
  spe,
  score_var = "SNAPa_all_genes_score",
  title = "SNAP-a\n(all genes)"
)

ggsave(
  filename = file.path(fld_plot, "spot_plot_SNAPa_all_genes_rep_samples.pdf"),
  plot = combined_plot_snapa_all_genes,
  height = 5.5, width = 11.5,
  units = "in"
)

## SNAP-n all-genes spot plot ----
combined_plot_snapn_all_genes <- plot_snap_spot(
  spe,
  score_var = "SNAPn_all_genes_score",
  title = "SNAP-n\n(all genes)"
)

ggsave(
  filename = file.path(fld_plot, "spot_plot_SNAPn_all_genes_rep_samples.pdf"),
  plot = combined_plot_snapn_all_genes,
  height = 5.5, width = 11.5,
  units = "in"
)

print("Finished making SNAP-a / SNAP-n top-2000 and all-genes spot plots for representative samples")


# Session Info ----
sessioninfo::session_info()
