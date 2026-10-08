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
fld_rds <- here("processed-data/rds/17_gsMap")
fld_plot <- here("plots/17_gsMap")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)

# Representative samples: Br8667 (NTC) and Br5973 (SCZ)
# (same pair used throughout the paper, e.g. 03_visium_spatial_clustering)
rep_sample_id <- c("V13M06-342_D1", "V13M06-343_D1")

dx_colors <- c(ntc = "blue", scz = "red")
dx_labels <- c(ntc = "NTC", scz = "SCZ")


# Load data ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)
spd_levels <- levels(spe$spd_label)

spot_res_df <- readRDS(file.path(fld_rds, "gsMap_spot_res_df.rds"))
sample_spd_df <- readRDS(file.path(fld_rds, "gsMap_sample_spd_df.rds")) |>
  mutate(
    spd_label = factor(spd_label, levels = spd_levels),
    dx = factor(dx, levels = c("ntc", "scz"))
  )
across_cauchy_df <- readRDS(file.path(fld_rds, "gsMap_across_sample_cauchy_df.rds")) |>
  mutate(spd_label = factor(spd_label, levels = spd_levels))

trait_names <- unique(spot_res_df$trait)

# Bonferroni threshold across spatial domains
bonf_line <- -log10(0.05 / length(spd_levels))


# Spatial domain color palette ----
spd_palette <- set_names(
  Polychrome::palette36.colors(length(spd_levels)),
  spd_levels
)


# 1) Spot plots of -log10(p) for representative samples ----------------------
## spot border = spatial domain, spot fill = gsMap spot-level -log10(p)
spe_rep <- spe[, spe$sample_id %in% rep_sample_id]
spe_rep$sample_label <- paste0(spe_rep$brnum, "_", toupper(spe_rep$dx))
rm(spe)

# error prevention
stopifnot(length(unique(spe_rep$sample_id)) == 2)

plot_gsmap_spot <- function(spe, trait) {
  trait_df <- spot_res_df |> filter(trait == !!trait)
  spe$neg_log10_p <- trait_df$neg_log10_p[match(spe$key, trait_df$key)]

  # error prevention
  stopifnot(!anyNA(spe$neg_log10_p))

  plot_list <- unique(spe$sample_label) |>
    set_names() |>
    map(.f = function(.smp) {
      make_escheR(spe[, spe$sample_label == .smp]) |>
        add_fill("neg_log10_p", point_size = 2.1) |>
        add_ground("spd_label", point_size = 2.2) +
        labs(title = .smp) +
        theme(
          plot.title = element_text(size = 20, hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)
        )
    })

  ## Shared color scale (grayscale: white = min, black = max) ----
  fill_scale <- scale_fill_gradient(
    name = paste0(trait, "\n-log10(p)"),
    limits = range(spe$neg_log10_p),
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

  ggpubr::ggarrange(
    plotlist = plot_list,
    nrow = 1, ncol = 2,
    common.legend = TRUE, legend = "right"
  )
}

for (.trait in trait_names) {
  ggsave(
    filename = file.path(fld_plot, paste0("spot_plot_", .trait, "_rep_samples.pdf")),
    plot = plot_gsmap_spot(spe_rep, .trait),
    height = 5.5, width = 11.5, units = "in"
  )
}


# 2) Across-sample Cauchy p per spatial domain (all / NTC / SCZ) ------------
p_across <- across_cauchy_df |>
  mutate(group = factor(group, levels = c("all", "ntc", "scz"))) |>
  ggplot(aes(x = spd_label, y = neg_log10_p_cauchy, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = bonf_line, linetype = "dashed") +
  scale_fill_manual(
    values = c(all = "grey40", dx_colors),
    labels = c(all = "All", dx_labels),
    name = NULL
  ) +
  facet_wrap(~trait, scales = "free_y") +
  labs(
    x = "Spatial domain",
    y = "-log10(Cauchy p), across samples",
    caption = "Dashed line: Bonferroni p = 0.05 across spatial domains"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(fld_plot, "across_sample_cauchy_by_spatial_domain.pdf"),
  p_across,
  height = 4.5, width = 4 + 3 * length(trait_names), units = "in"
)


# 3) Sample x domain heatmap of per-sample Cauchy p, split by diagnosis -----
for (.trait in trait_names) {
  heat_df <- sample_spd_df |>
    filter(trait == .trait) |>
    mutate(brnum = fct_reorder(brnum, neg_log10_p_cauchy, .fun = max))

  p_heat <- heat_df |>
    ggplot(aes(x = spd_label, y = brnum, fill = neg_log10_p_cauchy)) +
    geom_tile() +
    scale_fill_gradient(
      low = "white", high = "black",
      name = "-log10(Cauchy p)"
    ) +
    facet_wrap(~dx, scales = "free_y", labeller = as_labeller(dx_labels)) +
    labs(x = "Spatial domain", y = "Donor", title = .trait) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(fld_plot, paste0("sample_by_domain_cauchy_heatmap_", .trait, ".pdf")),
    p_heat,
    height = 9, width = 8, units = "in"
  )
}


# 4) Donor-level comparison of SCZ vs NTC per spatial domain -----------------
## Response: per-donor mean spot-level -log10(p) within a domain (less
## sensitive to spot counts than the per-sample Cauchy p).
## Model per trait x domain: mean_neg_log10_p ~ dx + age + sex
dx_test_df <- sample_spd_df |>
  group_by(trait, spd_label) |>
  group_modify(.f = function(.df, .key) {
    fit <- lm(mean_neg_log10_p ~ dx + age + sex, data = .df)
    coef(summary(fit))["dxscz", , drop = FALSE] |>
      as_tibble() |>
      set_names(c("estimate", "std_error", "t_value", "p_value")) |>
      mutate(n_donors = nrow(.df))
  }) |>
  group_by(trait) |>
  mutate(fdr = p.adjust(p_value, method = "BH")) |>
  ungroup()

write_csv(dx_test_df, file.path(fld_rds, "gsMap_dx_test_by_spatial_domain.csv"))
print(dx_test_df, n = Inf)

for (.trait in trait_names) {
  p_box <- sample_spd_df |>
    filter(trait == .trait) |>
    ggplot(aes(x = spd_label, y = mean_neg_log10_p, color = dx)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1) +
    scale_color_manual(values = dx_colors, labels = dx_labels, name = NULL) +
    labs(
      x = "Spatial domain",
      y = "Donor mean spot -log10(p)",
      title = .trait
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(fld_plot, paste0("donor_mean_neg_log10_p_by_dx_", .trait, ".pdf")),
    p_box,
    height = 4.5, width = 8, units = "in"
  )
}

print("Finished gsMap result visualizations")


# Session Info ----
sessioninfo::session_info()
