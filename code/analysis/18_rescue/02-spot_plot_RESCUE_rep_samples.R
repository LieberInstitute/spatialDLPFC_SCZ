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

source(here("code/analysis/18_rescue/rescue_plot_utils.R"))


# Config ----
fld_plot <- here("plots/18_rescue")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)


# Load Data ----
spe <- readRDS(
  here(
    "processed-data/rds/18_rescue",
    "fnl_spe_kept_spots_only_w_RESCUE.rds"
  )
)

ct_cols <- metadata(spe)$RESCUE_celltype_cols
stopifnot(length(ct_cols) > 0, all(ct_cols %in% colnames(colData(spe))))


# Subset to representative samples ----
# error prevention
stopifnot(all(rep_sample_id %in% unique(spe$sample_id)))

spe <- spe[, spe$sample_id %in% rep_sample_id]

# error prevention
stopifnot(length(unique(spe$sample_id)) == 2)
stopifnot("factor" %in% class(spe$spd_label))

## sample_label used for panel titles, e.g. Br8667_NTC ----
## (factor ordered as in `rep_sample_id`, i.e. NTC first)
spe$sample_label <- paste0(spe$brnum, "_", toupper(spe$dx))
spe$sample_label <- factor(
  spe$sample_label,
  levels = unique(spe$sample_label[match(rep_sample_id, spe$sample_id)])
)

## Broad-class proportions = sum of fine cell type proportions ----
broad_cols <- split(ct_cols, rescue_broad_class(ct_cols))
stopifnot(setequal(names(broad_cols), broad_class_levels))
for (.cls in names(broad_cols)) {
  spe[[paste0("RESCUE_broad_", .cls)]] <- rowSums(
    as.matrix(colData(spe)[, broad_cols[[.cls]], drop = FALSE])
  )
}
broad_vars <- paste0("RESCUE_broad_", broad_class_levels)


# Spatial domain color palette ----
spd_palette <- set_names(
  Polychrome::palette36.colors(length(levels(spe$spd_label))),
  levels(spe$spd_label)
)


# Spot plots -----------------------------------------------------------------
## escheR spot plot of one proportion column for both samples, sharing a
## [0, max] grayscale fill (white = 0, black = max across both samples) ----
plot_prop_spot <- function(spe, prop_var, title) {
  fill_scale <- scale_fill_gradient(
    name = "Proportion",
    limits = c(0, max(spe[[prop_var]], na.rm = TRUE)),
    low = "white", high = "black"
  )

  plot_list <- levels(spe$sample_label) |>
    set_names() |>
    map(.f = function(.smp) {
      sub_spe <- spe[, spe$sample_label == .smp]

      make_escheR(sub_spe) |>
        add_fill(prop_var, point_size = 2.1) |>
        add_ground("spd_label", point_size = 2.2) +
        fill_scale +
        scale_color_manual(
          name = "Spatial Domain",
          values = spd_palette,
          guide = guide_legend(override.aes = list(size = 7))
        ) +
        labs(title = .smp) +
        theme(
          plot.title = element_text(size = 20, hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)
        )
    })

  ggpubr::ggarrange(
    plotlist = plot_list,
    nrow = 1, ncol = 2,
    common.legend = TRUE,
    legend = "right"
  ) |>
    ggpubr::annotate_figure(
      top = ggpubr::text_grob(title, size = 22, face = "bold")
    )
}

## Broad classes: one page per class ----
pdf(
  file.path(fld_plot, "spot_plot_RESCUE_broad_class_rep_samples.pdf"),
  height = 6, width = 11.5,
  bg = "white"
)
walk2(broad_vars, broad_class_levels, \(.var, .cls) {
  print(plot_prop_spot(spe, .var, sprintf("RESCUE: %s (broad class)", .cls)))
})
dev.off()

## Fine cell types: one page per cell type ----
pdf(
  file.path(fld_plot, "spot_plot_RESCUE_celltype_rep_samples.pdf"),
  height = 6, width = 11.5,
  bg = "white"
)
walk(ct_cols, \(.var) {
  print(plot_prop_spot(spe, .var, sprintf("RESCUE: %s", str_remove(.var, "^RESCUE_"))))
})
dev.off()


# Composition per spatial domain within each representative sample ---------
rep_df <- colData(spe) |>
  as.data.frame() |>
  select(key, sample_label, spd_label, all_of(ct_cols))

## Stacked bar: mean broad-class composition per spatial domain ----
broad_comp_df <- rep_df |>
  pivot_longer(all_of(ct_cols), names_to = "cell_type", values_to = "prop") |>
  mutate(
    broad_class = factor(rescue_broad_class(cell_type), levels = broad_class_levels)
  ) |>
  # sum fine types within spot, then average across spots in a domain
  group_by(sample_label, spd_label, key, broad_class) |>
  summarize(prop = sum(prop), .groups = "drop") |>
  group_by(sample_label, spd_label, broad_class) |>
  summarize(mean_prop = mean(prop), .groups = "drop")

# error prevention: mean compositions should still sum to 1
stopifnot(
  broad_comp_df |>
    group_by(sample_label, spd_label) |>
    summarize(s = sum(mean_prop), .groups = "drop") |>
    pull(s) |>
    (\(s) all(abs(s - 1) < 1e-6))()
)

## Both versions: incl. unattributed, and conditional (excl. unattributed)
broad_comp_list <- list(
  with_unattributed = broad_comp_df,
  excl_unattributed = to_conditional_comp(
    broad_comp_df, c("sample_label", "spd_label")
  )
)

pwalk(comp_versions, \(version, file_suffix, y_label) {
  broad_bar_plot <- ggplot(
    broad_comp_list[[version]],
    aes(x = spd_label, y = mean_prop, fill = fct_rev(broad_class))
  ) +
    geom_col(width = 0.8, color = "white", linewidth = 0.3) +
    facet_wrap(~sample_label, nrow = 1) +
    scale_fill_broad_class() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Spatial Domain", y = y_label) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12, face = "bold")
    )

  ggsave(
    file.path(
      fld_plot,
      sprintf("RESCUE_broad_class_by_spd_rep_samples%s.pdf", file_suffix)
    ),
    broad_bar_plot,
    width = 9, height = 5, units = "in"
  )
})

## Heatmap: mean fine cell type proportion per spatial domain ----
fine_comp_df <- rep_df |>
  pivot_longer(all_of(ct_cols), names_to = "cell_type", values_to = "prop") |>
  group_by(sample_label, spd_label, cell_type) |>
  summarize(mean_prop = mean(prop), .groups = "drop") |>
  mutate(
    cell_type = factor(str_remove(cell_type, "^RESCUE_"), levels = rev(str_remove(ct_cols, "^RESCUE_")))
  )

fine_heatmap_plot <- ggplot(
  fine_comp_df,
  aes(x = spd_label, y = cell_type, fill = mean_prop)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~sample_label, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#0d366b", name = "Mean\nproportion") +
  labs(x = "Spatial Domain", y = "RESCUE cell type") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(fld_plot, "RESCUE_celltype_by_spd_heatmap_rep_samples.pdf"),
  fine_heatmap_plot,
  width = 8, height = 8.5, units = "in"
)

print("Finished making RESCUE spot plots and per-domain compositions for representative samples")


# Session Info ----
sessioninfo::session_info()
