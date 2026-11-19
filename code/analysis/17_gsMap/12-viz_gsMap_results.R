# Load Packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  library(here)
  library(escheR)
  library(Polychrome)
  library(ggpubr)
  library(ggridges)
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
  sort(spd_levels) # assign by SpD number, not laminar order
)


# 1) Spot plots of -log10(p) for representative samples ----------------------
## spot border = spatial domain, spot fill = gsMap spot-level -log10(p)
spe_rep <- spe[, spe$sample_id %in% rep_sample_id]
spe_rep$sample_label <- paste0(spe_rep$brnum, "_", toupper(spe_rep$dx))

# spot-level metadata for the per-domain ridge plots (section 5)
spot_meta <- colData(spe) |>
  as.data.frame() |>
  select(key, brnum, dx, spd_label) |>
  as_tibble()
rm(spe)

# error prevention
stopifnot(length(unique(spe_rep$sample_id)) == 2)

plot_gsmap_spot <- function(spe, trait) {
  trait_df <- spot_res_df |> filter(trait == !!trait)
  spe$neg_log10_p <- trait_df$neg_log10_p[match(spe$key, trait_df$key)]

  # Spots in sample x domain combinations with < 30 spots have no gsMap
  # result (dropped by gsMap; see 11-collect_gsMap_results.R) and are shown
  # with the NA colour below.
  n_na <- sum(is.na(spe$neg_log10_p))
  if (n_na > 0) {
    cat(sprintf("%s: %d spots without gsMap result (small domains)\n", trait, n_na))
  }

  smp_labels <- unique(spe$sample_label) |> set_names()

  spot_theme <- theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )

  ## Shared grayscale (white = min, black = max) across the two samples
  fill_scale <- scale_fill_gradient(
    name = paste0(trait, "\n-log10(p)"),
    limits = range(spe$neg_log10_p, na.rm = TRUE),
    low = "white", high = "black",
    # grey50 (the default) would blend into the grayscale; use a distinct colour
    na.value = "lightskyblue"
  )

  ## Row 1: composite, spot fill = -log10(p), spot border = spatial domain ----
  combo_list <- smp_labels |>
    map(.f = function(.smp) {
      make_escheR(spe[, spe$sample_label == .smp]) |>
        add_fill("neg_log10_p", point_size = 2.1) |>
        add_ground("spd_label", point_size = 2.2) +
        fill_scale +
        scale_color_manual(name = "Spatial Domain", values = spd_palette) +
        labs(title = .smp) +
        spot_theme
    })

  ## Row 2: spatial domains only ----
  spd_list <- smp_labels |>
    map(.f = function(.smp) {
      make_escheR(spe[, spe$sample_label == .smp]) |>
        add_fill("spd_label", point_size = 2.1) +
        scale_fill_manual(
          name = "Spatial Domain",
          values = spd_palette,
          guide = guide_legend(override.aes = list(size = 7))
        ) +
        labs(title = .smp) +
        spot_theme
    })

  ## Row 3: gsMap -log10(p) only ----
  p_list <- smp_labels |>
    map(.f = function(.smp) {
      make_escheR(spe[, spe$sample_label == .smp]) |>
        add_fill("neg_log10_p", point_size = 2.1) +
        fill_scale +
        labs(title = .smp) +
        spot_theme
    })

  ## Shared legends at the bottom: spatial domain (row 2) and -log10(p) (row 3)
  legend_bottom <- function(.p) {
    ggpubr::get_legend(.p + theme(legend.position = "bottom"))
  }
  legend_row <- ggpubr::ggarrange(
    legend_bottom(
      spd_list[[1]] +
        guides(fill = guide_legend(nrow = 2, override.aes = list(size = 5)))
    ),
    legend_bottom(
      p_list[[1]] +
        guides(fill = guide_colorbar(barwidth = 15, title.vjust = 0.8))
    ),
    nrow = 2, ncol = 1
  )

  ggpubr::ggarrange(
    plotlist = c(combo_list, spd_list, p_list),
    nrow = 3, ncol = 2,
    legend = "none"
  ) |>
    ggpubr::ggarrange(
      legend_row,
      nrow = 2, ncol = 1, heights = c(15, 2)
    )
}

for (.trait in trait_names) {
  ggsave(
    filename = file.path(fld_plot, paste0("spot_plot_", .trait, "_rep_samples.pdf")),
    plot = plot_gsmap_spot(spe_rep, .trait),
    height = 17, width = 11, units = "in"
  )
}


# 2) Across-sample Cauchy p per spatial domain (all / NTC / SCZ) ------------
## One facet per group (all / NTC / SCZ); rows = traits
group_labels <- c(
  all = "All (n = 63)",
  ntc = "NTC (n = 31)",
  scz = "SCZ (n = 32)"
)

p_across <- across_cauchy_df |>
  mutate(group = factor(group, levels = c("all", "ntc", "scz"))) |>
  ggplot(aes(x = spd_label, y = neg_log10_p_cauchy, fill = group)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_hline(yintercept = bonf_line, linetype = "dashed") +
  scale_fill_manual(values = c(all = "grey40", dx_colors)) +
  facet_grid(trait ~ group, labeller = labeller(group = group_labels)) +
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
  height = 1 + 3.5 * length(trait_names), width = 10, units = "in"
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

  ## Dot-plot version: size = number of spots, colour = -log10(Cauchy p) ----
  ## Same rows (donor order), columns and facets as the heatmap. Donor x
  ## domain pairs without a gsMap result (domain absent, or < 30 spots and
  ## dropped by gsMap) are marked with a cross.
  dot_df <- heat_df |>
    distinct(sample_id, brnum, dx) |>
    expand_grid(spd_label = factor(spd_levels, levels = spd_levels)) |>
    left_join(
      heat_df |> select(sample_id, spd_label, n_spots, neg_log10_p_cauchy),
      by = c("sample_id", "spd_label")
    ) |>
    mutate(missing = is.na(neg_log10_p_cauchy))

  # error prevention: one row per donor x domain
  stopifnot(nrow(dot_df) == n_distinct(heat_df$sample_id) * length(spd_levels))
  cat(sprintf(
    "%s: %d donor x domain pairs without gsMap result\n",
    .trait, sum(dot_df$missing)
  ))

  p_dot <- ggplot(dot_df, aes(x = spd_label, y = brnum)) +
    geom_point(
      data = \(.d) filter(.d, !missing),
      aes(size = n_spots, fill = neg_log10_p_cauchy),
      shape = 21, colour = "grey30", stroke = 0.2
    ) +
    geom_point(
      data = \(.d) filter(.d, missing),
      aes(shape = "No gsMap result\n(domain absent or < 30 spots)"),
      size = 2.5, colour = "red3", stroke = 0.8
    ) +
    scale_fill_gradient(
      low = "white", high = "black",
      name = "-log10(Cauchy p)"
    ) +
    scale_size_area(name = "Number of spots", max_size = 6) +
    scale_shape_manual(values = 4, name = NULL) +
    guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(order = 2),
      shape = guide_legend(order = 3)
    ) +
    facet_wrap(~dx, scales = "free_y", labeller = as_labeller(dx_labels)) +
    labs(x = "Spatial domain", y = "Donor", title = .trait) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(fld_plot, paste0("sample_by_domain_cauchy_dotplot_", .trait, ".pdf")),
    p_dot,
    height = 9, width = 8.5, units = "in"
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
  box_df <- sample_spd_df |> filter(trait == .trait)

  ## Number of donors per diagnosis x domain (fewer where gsMap dropped a
  ## domain with < 30 spots or the domain is absent; see section 3)
  n_df <- box_df |>
    count(dx, spd_label, name = "n_donors") |>
    mutate(y = max(box_df$mean_neg_log10_p) * 1.08)

  ## Facet strips show the total number of donors per diagnosis
  facet_labels <- box_df |>
    distinct(dx, sample_id) |>
    count(dx) |>
    mutate(label = sprintf("%s (n = %d donors)", dx_labels[as.character(dx)], n)) |>
    select(dx, label) |>
    deframe()

  set.seed(20241118) # reproducible jitter
  p_box <- box_df |>
    ggplot(aes(x = spd_label, y = mean_neg_log10_p, color = dx)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.15, height = 0, size = 1) +
    geom_text(
      data = n_df,
      aes(y = y, label = paste0("n=", n_donors)),
      size = 3, show.legend = FALSE
    ) +
    scale_color_manual(values = dx_colors, guide = "none") +
    facet_wrap(~dx, nrow = 1, labeller = as_labeller(facet_labels)) +
    labs(
      x = "Spatial domain",
      y = "Donor mean spot -log10(p)",
      title = .trait,
      caption = "n = number of donors with a gsMap result in the domain"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(fld_plot, paste0("donor_mean_neg_log10_p_by_dx_", .trait, ".pdf")),
    p_box,
    height = 4.5, width = 10, units = "in"
  )
}


# 5) Ridge plots of spot-level -log10(p), one page per spatial domain --------
## One ridge per donor (density over the donor's spots in that domain),
## filled by diagnosis; ridge label = donor (number of spots), coloured by
## diagnosis. Donors are ordered by their number of spots in the domain
## (most spots on top).
## Donor x domain pairs without gsMap results (see section 3) have no ridge.
ridge_df <- spot_res_df |>
  select(key, trait, neg_log10_p) |>
  inner_join(spot_meta, by = "key", relationship = "many-to-one") |>
  mutate(dx = factor(dx, levels = c("ntc", "scz")))

# error prevention: every spot result has metadata
stopifnot(nrow(ridge_df) == nrow(spot_res_df))

for (.trait in trait_names) {
  pdf(
    file.path(fld_plot, paste0("ridge_neg_log10_p_by_spatial_domain_", .trait, ".pdf")),
    height = 12, width = 7
  )

  for (.spd in spd_levels) {
    .df <- ridge_df |>
      filter(trait == .trait, spd_label == .spd) |>
      group_by(brnum) |>
      mutate(
        n_spots = n(),
        ridge_label = sprintf("%s (n = %s)", brnum, format(n_spots, big.mark = ","))
      ) |>
      ungroup() |>
      # most spots on top (first factor level is drawn at the bottom)
      arrange(n_spots, brnum) |>
      mutate(ridge_label = fct_inorder(ridge_label))

    # y-axis label colours by diagnosis, in factor-level order
    label_colors <- .df |>
      distinct(ridge_label, dx) |>
      arrange(ridge_label) |>
      pull(dx) |>
      as.character()
    label_colors <- dx_colors[label_colors]
    stopifnot(length(label_colors) == nlevels(.df$ridge_label))

    p_ridge <- ggplot(.df, aes(x = neg_log10_p, y = ridge_label, fill = dx)) +
      geom_density_ridges(
        scale = 2.5, alpha = 0.6, linewidth = 0.25,
        rel_min_height = 0.005, from = 0
      ) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
      scale_fill_manual(values = dx_colors, labels = dx_labels, name = NULL) +
      scale_y_discrete(expand = expansion(add = c(0.3, 2.5))) +
      labs(
        x = "Spot -log10(p)",
        y = "Donor (number of spots)",
        title = paste0(.trait, ": ", .spd),
        subtitle = sprintf(
          "%s spots from %d donors",
          format(nrow(.df), big.mark = ","), n_distinct(.df$brnum)
        ),
        caption = "Dashed line: p = 0.05"
      ) +
      theme_ridges(font_size = 9, grid = TRUE) +
      theme(
        legend.position = "top",
        # vectorised colours work but are not officially supported by ggplot2
        axis.text.y = suppressWarnings(element_text(colour = label_colors))
      )

    print(p_ridge)
  }

  dev.off()
}


# 6) Per-donor spot plots: spatial domain | -log10(p), one page per donor ----
## The -log10(p) colour scale is shared across all donors (per trait), so
## pages are directly comparable. Spots in domains dropped by gsMap (< 30
## spots) have no result and are shown in light blue.
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

# NTC first, then SCZ; donors in BrNum order within each group
donor_df <- colData(spe) |>
  as.data.frame() |>
  distinct(sample_id, brnum, dx) |>
  mutate(dx = factor(dx, levels = c("ntc", "scz"))) |>
  arrange(dx, brnum)

# error prevention: one section per donor
stopifnot(!anyDuplicated(donor_df$brnum), !anyDuplicated(donor_df$sample_id))

donor_theme <- theme(
  plot.title = element_text(size = 14, hjust = 0.5),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
)

for (.trait in trait_names) {
  trait_df <- spot_res_df |> filter(trait == .trait)
  spe$neg_log10_p <- trait_df$neg_log10_p[match(spe$key, trait_df$key)]

  fill_scale <- scale_fill_gradient(
    name = paste0(.trait, "\n-log10(p)"),
    limits = range(spe$neg_log10_p, na.rm = TRUE),
    low = "white", high = "black",
    na.value = "lightskyblue"
  )

  pdf(
    file.path(fld_plot, paste0("spot_plot_", .trait, "_all_donors.pdf")),
    height = 6, width = 13
  )

  for (.i in seq_len(nrow(donor_df))) {
    .smp <- donor_df$sample_id[.i]
    .spe <- spe[, spe$sample_id == .smp]
    n_na <- sum(is.na(.spe$neg_log10_p))

    p_spd <- make_escheR(.spe) |>
      add_fill("spd_label", point_size = 1.6) +
      scale_fill_manual(
        name = "Spatial Domain",
        values = spd_palette,
        drop = FALSE,
        guide = guide_legend(override.aes = list(size = 5))
      ) +
      labs(title = "Spatial domain") +
      donor_theme

    p_nlp <- make_escheR(.spe) |>
      add_fill("neg_log10_p", point_size = 1.6) +
      fill_scale +
      labs(title = "Spot -log10(p)") +
      donor_theme

    p_page <- ggpubr::ggarrange(p_spd, p_nlp, nrow = 1, ncol = 2) |>
      ggpubr::annotate_figure(
        top = ggpubr::text_grob(
          sprintf(
            "%s (%s) | %s | %s spots%s",
            donor_df$brnum[.i], dx_labels[as.character(donor_df$dx[.i])],
            .smp, format(ncol(.spe), big.mark = ","),
            if (n_na > 0) sprintf(", %d without gsMap result (light blue)", n_na) else ""
          ),
          size = 16, face = "bold"
        )
      )

    suppressMessages(print(p_page))
  }

  dev.off()
}

rm(spe)

print("Finished gsMap result visualizations")


# Session Info ----
sessioninfo::session_info()
