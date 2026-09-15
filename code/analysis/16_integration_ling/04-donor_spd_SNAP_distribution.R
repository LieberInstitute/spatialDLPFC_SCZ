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
    cols = c(SNAPa_zscore, SNAPn_zscore),
    names_to = "program",
    values_to = "zscore"
  ) |>
  mutate(
    program = recode(
      program,
      "SNAPa_zscore" = "SNAP-a",
      "SNAPn_zscore" = "SNAP-n"
    ),
    dx = factor(dx, levels = c("ntc", "scz"))
  )

dx_colors <- c(ntc = "blue", scz = "red")
dx_labels <- c(ntc = "NTC", scz = "SCZ")


# 1) Distribution across donors ---------------------------------------------
## Order donors within each dx group by median normalized score, pooled
## across both programs so the ordering is shared across facets ----
donor_summary <- score_long |>
  group_by(brnum, dx) |>
  summarize(med_score = median(zscore), .groups = "drop") |>
  arrange(dx, med_score)

# error prevention: one ordering row per donor
stopifnot(!anyDuplicated(donor_summary$brnum))

score_long <- score_long |>
  mutate(brnum = factor(brnum, levels = donor_summary$brnum))

donor_plot <- ggplot(
  score_long,
  aes(x = brnum, y = zscore, fill = dx)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~program, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = dx_colors, labels = dx_labels, name = "Diagnosis") +
  labs(
    x = "Donor",
    y = "Normalized SNAP projection score (z)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(fld_plot, "SNAP_score_by_donor.pdf"),
  donor_plot,
  width = 13, height = 7, units = "in"
)


# 2) Distribution across spatial domains, split by diagnosis -----------------
domain_plot <- ggplot(
  score_long,
  aes(x = spd_label, y = zscore, fill = dx)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
  facet_wrap(~program, ncol = 1) +
  scale_fill_manual(values = dx_colors, labels = dx_labels, name = "Diagnosis") +
  labs(
    x = "Spatial Domain",
    y = "Normalized SNAP projection score (z)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(fld_plot, "SNAP_score_by_spatial_domain.pdf"),
  domain_plot,
  width = 9, height = 7, units = "in"
)


# 3) Donor x spatial-domain mean-score heatmap -------------------------------
donor_domain_mean <- score_long |>
  group_by(program, brnum, dx, spd_label) |>
  summarize(mean_zscore = mean(zscore), .groups = "drop")

heatmap_plot <- ggplot(
  donor_domain_mean,
  aes(x = spd_label, y = brnum, fill = mean_zscore)
) +
  geom_tile(color = "grey90") +
  facet_grid(dx ~ program, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, name = "Mean\nz-score"
  ) +
  labs(x = "Spatial Domain", y = "Donor") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave(
  file.path(fld_plot, "SNAP_score_donor_by_domain_heatmap.pdf"),
  heatmap_plot,
  width = 8.5, height = 11, units = "in"
)

print("Finished making donor-level SNAP-a / SNAP-n distribution plots")


# Session info ----
sessioninfo::session_info()
