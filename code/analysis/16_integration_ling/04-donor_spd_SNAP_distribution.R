# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})


# Config ----
fld_plot <- here("plots/16_integration_ling")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)


# Load data ----
score_df <- readRDS(
  here(
    "processed-data/rds/16_integration_ling",
    "SNAP_score_df.rds"
  )
)


# Reshape to long format ----
score_long <- score_df |>
  pivot_longer(
    cols = c(
      SNAPa_top_2000_score, SNAPn_top_2000_score,
      SNAPa_all_genes_score, SNAPn_all_genes_score
    ),
    names_to = "score_name",
    values_to = "score"
  ) |>
  mutate(
    program = if_else(str_starts(score_name, "SNAPa"), "SNAP-a", "SNAP-n"),
    gene_set = if_else(str_detect(score_name, "top_2000"), "top_2000", "all_genes"),
    dx = factor(dx, levels = c("ntc", "scz"))
  ) |>
  select(-score_name)

dx_colors <- c(ntc = "blue", scz = "red")
dx_labels <- c(ntc = "NTC", scz = "SCZ")

# Gene sets are plotted separately below, each producing its own set of
# donor-, spatial-domain-, and heatmap-level plots.
gene_set_info <- tribble(
  ~gene_set, ~label, ~file_suffix,
  "top_2000", "top 2000", "top_2000",
  "all_genes", "all genes", "all_genes"
)


# Make donor / spatial-domain / heatmap plots for one gene set --------------
make_snap_distribution_plots <- function(score_long, gene_set_name, label, file_suffix) {
  score_sub <- score_long |> filter(gene_set == gene_set_name)
  y_label <- sprintf("SNAP projection score (%s)", label)

  ## 1) Distribution across donors ----
  ## Order donors within each dx group by median score, pooled across both
  ## programs so the ordering is shared across facets ----
  donor_summary <- score_sub |>
    group_by(brnum, dx) |>
    summarize(med_score = median(score), .groups = "drop") |>
    arrange(dx, med_score)

  # error prevention: one ordering row per donor
  stopifnot(!anyDuplicated(donor_summary$brnum))

  score_sub <- score_sub |>
    mutate(brnum = factor(brnum, levels = donor_summary$brnum))

  donor_plot <- ggplot(
    score_sub,
    aes(x = brnum, y = score, fill = dx)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    facet_wrap(~program, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = dx_colors, labels = dx_labels, name = "Diagnosis") +
    labs(
      x = "Donor",
      y = y_label
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      strip.text = element_text(size = 12, face = "bold")
    )

  ggsave(
    file.path(fld_plot, sprintf("SNAP_score_by_donor_%s.pdf", file_suffix)),
    donor_plot,
    width = 13, height = 7, units = "in"
  )

  ## 2) Distribution across spatial domains, split by diagnosis ----
  domain_plot <- ggplot(
    score_sub,
    aes(x = spd_label, y = score, fill = dx)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
    facet_wrap(~program, ncol = 1) +
    scale_fill_manual(values = dx_colors, labels = dx_labels, name = "Diagnosis") +
    labs(
      x = "Spatial Domain",
      y = y_label
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12, face = "bold")
    )

  ggsave(
    file.path(fld_plot, sprintf("SNAP_score_by_spatial_domain_%s.pdf", file_suffix)),
    domain_plot,
    width = 9, height = 7, units = "in"
  )

  ## 3) Donor x spatial-domain mean-score heatmap ----
  donor_domain_mean <- score_sub |>
    group_by(program, brnum, dx, spd_label) |>
    summarize(mean_score = mean(score), .groups = "drop")

  heatmap_plot <- ggplot(
    donor_domain_mean,
    aes(x = spd_label, y = brnum, fill = mean_score)
  ) +
    geom_tile(color = "grey90") +
    facet_grid(dx ~ program, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = 0, name = "Mean\nscore"
    ) +
    labs(x = "Spatial Domain", y = "Donor") +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 11, face = "bold")
    )

  ggsave(
    file.path(fld_plot, sprintf("SNAP_score_donor_by_domain_heatmap_%s.pdf", file_suffix)),
    heatmap_plot,
    width = 8.5, height = 11, units = "in"
  )
}

pwalk(
  gene_set_info,
  \(gene_set, label, file_suffix) {
    make_snap_distribution_plots(score_long, gene_set, label, file_suffix)
  }
)

print("Finished making donor-level SNAP-a / SNAP-n distribution plots (top-2000 and all-genes)")


# Session info ----
sessioninfo::session_info()
