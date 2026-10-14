# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

source(here("code/analysis/18_rescue/rescue_plot_utils.R"))


# Config ----
fld_plot <- here("plots/18_rescue")
dir.create(fld_plot, showWarnings = FALSE, recursive = TRUE)

# Donor x domain combinations with fewer spots than this are excluded from
# plots (their mean proportion is too noisy); they are kept in the saved csv.
min_spots <- 10


# Load data ----
## Slim key-indexed RESCUE proportions, from 01-merge_RESCUE_to_spe.R ----
prop_df <- readRDS(
  here("processed-data/rds/18_rescue", "RESCUE_prop_df.rds")
)

ct_cols <- setdiff(
  grep("^RESCUE_", colnames(prop_df), value = TRUE),
  c("RESCUE_nUMI", "RESCUE_unattributed_ratio")
)

# error prevention: spots without a RESCUE result are dropped
n_na <- sum(is.na(prop_df$RESCUE_Astro))
cat(sprintf("Dropping %d spots without RESCUE proportions\n", n_na))
prop_df <- prop_df |> filter(!is.na(RESCUE_Astro))


# Donor x spatial domain mean proportions ------------------------------------
## Each donor contributes one sample; spot proportions are averaged within
## donor x domain (and donor overall) so downstream summaries weight donors,
## not spots, equally.
prop_long <- prop_df |>
  pivot_longer(all_of(ct_cols), names_to = "cell_type", values_to = "prop") |>
  mutate(
    cell_type = str_remove(cell_type, "^RESCUE_"),
    broad_class = factor(rescue_broad_class(cell_type), levels = broad_class_levels),
    cell_type = factor(cell_type, levels = str_remove(ct_cols, "^RESCUE_")),
    dx = factor(dx, levels = c("ntc", "scz"))
  )

donor_spd_fine <- prop_long |>
  group_by(brnum, dx, spd_label, cell_type, broad_class) |>
  summarize(mean_prop = mean(prop), n_spots = n(), .groups = "drop")

donor_all_fine <- prop_long |>
  group_by(brnum, dx, cell_type, broad_class) |>
  summarize(mean_prop = mean(prop), .groups = "drop") |>
  mutate(spd_label = "All spots")

donor_spd_broad <- bind_rows(donor_spd_fine, donor_all_fine) |>
  mutate(
    spd_label = factor(spd_label, levels = c(levels(prop_df$spd_label), "All spots"))
  ) |>
  group_by(brnum, dx, spd_label, broad_class) |>
  summarize(mean_prop = sum(mean_prop), .groups = "drop")

# error prevention: broad-class means sum to 1 within each donor x domain
stopifnot(
  donor_spd_broad |>
    group_by(brnum, spd_label) |>
    summarize(s = sum(mean_prop), .groups = "drop") |>
    pull(s) |>
    (\(s) all(abs(s - 1) < 1e-6))()
)

write_csv(
  donor_spd_fine,
  here("processed-data/rds/18_rescue", "RESCUE_donor_spd_mean_prop.csv")
)

## Drop sparse donor x domain combinations from all plots below ----
sparse_df <- donor_spd_fine |>
  distinct(brnum, dx, spd_label, n_spots) |>
  filter(n_spots < min_spots)
cat(sprintf(
  "Excluding %d donor x domain combinations with < %d spots from plots:\n",
  nrow(sparse_df), min_spots
))
print(sparse_df, n = Inf)

donor_spd_fine <- donor_spd_fine |>
  anti_join(sparse_df, by = c("brnum", "spd_label"))
donor_spd_broad <- donor_spd_broad |>
  anti_join(sparse_df, by = c("brnum", "spd_label"))


# Compositional (stacked bar) plots -------------------------------------------
## Each plot is made in two versions (see `comp_versions`): proportions as
## reported by RESCUE (incl. unattributed), and the conditional composition
## of attributed signal (excl. unattributed; renormalized per donor x domain
## before averaging across donors).
donor_spd_broad_list <- list(
  with_unattributed = donor_spd_broad,
  excl_unattributed = to_conditional_comp(
    donor_spd_broad, c("brnum", "spd_label")
  )
)

# error prevention: conditional compositions sum to 1 within donor x domain
stopifnot(
  donor_spd_broad_list$excl_unattributed |>
    group_by(brnum, spd_label) |>
    summarize(s = sum(mean_prop), .groups = "drop") |>
    pull(s) |>
    (\(s) all(abs(s - 1) < 1e-6))()
)

## Donors ordered within diagnosis by overall Excit proportion (incl.
## unattributed), shared by both versions of the donor plot ----
donor_order <- donor_spd_broad |>
  filter(spd_label == "All spots", broad_class == "Excit") |>
  arrange(dx, mean_prop) |>
  pull(brnum)

# error prevention: one ordering row per donor
stopifnot(!anyDuplicated(donor_order))

make_composition_plots <- function(version, file_suffix, y_label) {
  broad_df <- donor_spd_broad_list[[version]]

  ## 1) Overall composition per spatial domain (diagnoses pooled) ----
  ## Mean of donor means across all 63 donors ----
  spd_all_broad <- broad_df |>
    group_by(spd_label, broad_class) |>
    summarize(mean_prop = mean(mean_prop), .groups = "drop")

  spd_all_plot <- ggplot(
    spd_all_broad,
    aes(x = spd_label, y = mean_prop, fill = fct_rev(broad_class))
  ) +
    geom_col(width = 0.8, color = "white", linewidth = 0.3) +
    scale_fill_broad_class() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Spatial Domain", y = paste0(y_label, "\n(mean of donor means)")) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(fld_plot, sprintf("RESCUE_broad_class_by_spd%s.pdf", file_suffix)),
    spd_all_plot,
    width = 7, height = 5, units = "in"
  )

  ## 2) Composition per spatial domain, by diagnosis ----
  spd_dx_broad <- broad_df |>
    group_by(dx, spd_label, broad_class) |>
    summarize(mean_prop = mean(mean_prop), .groups = "drop")

  spd_dx_plot <- ggplot(
    spd_dx_broad,
    aes(x = dx, y = mean_prop, fill = fct_rev(broad_class))
  ) +
    geom_col(width = 0.8, color = "white", linewidth = 0.3) +
    facet_wrap(~spd_label, nrow = 1, strip.position = "bottom") +
    scale_x_discrete(labels = dx_labels) +
    scale_fill_broad_class() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = paste0(y_label, "\n(mean of donor means)")) +
    theme_classic(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.spacing = unit(4, "pt")
    )

  ggsave(
    file.path(fld_plot, sprintf("RESCUE_broad_class_by_spd_dx%s.pdf", file_suffix)),
    spd_dx_plot,
    width = 10, height = 5, units = "in"
  )

  ## 3) Composition per donor, stratified by spatial domain ----
  donor_plot <- broad_df |>
    mutate(brnum = factor(brnum, levels = donor_order)) |>
    ggplot(aes(x = brnum, y = mean_prop, fill = fct_rev(broad_class))) +
    geom_col(width = 0.9) +
    facet_grid(spd_label ~ dx, scales = "free_x", space = "free_x",
      labeller = labeller(dx = dx_labels)) +
    scale_fill_broad_class() +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    labs(x = "Donor", y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      strip.text.y = element_text(angle = 0, hjust = 0),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    file.path(fld_plot, sprintf("RESCUE_broad_class_by_donor_spd%s.pdf", file_suffix)),
    donor_plot,
    width = 14, height = 11, units = "in"
  )
}

pwalk(comp_versions, make_composition_plots)


# Fine cell type x spatial domain heatmap, by diagnosis --------------------
spd_dx_fine <- donor_spd_fine |>
  group_by(dx, spd_label, cell_type) |>
  summarize(mean_prop = mean(mean_prop), .groups = "drop") |>
  mutate(cell_type = fct_rev(cell_type))

fine_heatmap_plot <- ggplot(
  spd_dx_fine,
  aes(x = spd_label, y = cell_type, fill = mean_prop)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~dx, nrow = 1, labeller = labeller(dx = dx_labels)) +
  scale_fill_gradient(low = "white", high = "#0d366b", name = "Mean\nproportion") +
  labs(x = "Spatial Domain", y = "RESCUE cell type") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(fld_plot, "RESCUE_celltype_by_spd_heatmap_dx.pdf"),
  fine_heatmap_plot,
  width = 8, height = 8.5, units = "in"
)

## Row-scaled version: where in the cortex is each cell type enriched? ----
fine_heatmap_scaled_plot <- spd_dx_fine |>
  group_by(dx, cell_type) |>
  mutate(rel_prop = mean_prop / max(mean_prop)) |>
  ungroup() |>
  ggplot(aes(x = spd_label, y = cell_type, fill = rel_prop)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~dx, nrow = 1, labeller = labeller(dx = dx_labels)) +
  scale_fill_gradient(
    low = "white", high = "#0d366b", limits = c(0, 1),
    name = "Mean proportion\n(/ row max)"
  ) +
  labs(x = "Spatial Domain", y = "RESCUE cell type") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(fld_plot, "RESCUE_celltype_by_spd_heatmap_dx_row_scaled.pdf"),
  fine_heatmap_scaled_plot,
  width = 8, height = 8.5, units = "in"
)


# Donor-level distribution per spatial domain, by diagnosis ----------------
## One point per donor x domain; one page per broad class (fine types faceted)
plot_donor_box <- function(df, facet_var, title) {
  ggplot(df, aes(x = spd_label, y = mean_prop, fill = dx)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
    geom_point(
      aes(color = dx),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size = 0.6, alpha = 0.6, show.legend = FALSE
    ) +
    facet_wrap(vars({{ facet_var }}), scales = "free_y") +
    scale_fill_manual(values = dx_colors, labels = dx_labels, name = "Diagnosis") +
    scale_color_manual(values = dx_colors) +
    labs(x = "Spatial Domain", y = "Donor mean RESCUE proportion", title = title) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )
}

pdf(
  file.path(fld_plot, "RESCUE_donor_prop_by_spd_dx_boxplot.pdf"),
  width = 12, height = 8, bg = "white"
)

## Broad classes on the first page ----
print(plot_donor_box(
  donor_spd_broad |> filter(spd_label != "All spots"),
  broad_class,
  "RESCUE broad classes"
))

## Fine cell types, one page per broad class with >1 cell type ----
walk(broad_class_levels, \(.cls) {
  sub_df <- donor_spd_fine |> filter(broad_class == .cls)
  if (n_distinct(sub_df$cell_type) > 1) {
    print(plot_donor_box(sub_df, cell_type, sprintf("RESCUE cell types: %s", .cls)))
  }
})
dev.off()

print("Finished making donor / spatial domain stratified RESCUE composition plots")


# Session info ----
sessioninfo::session_info()
