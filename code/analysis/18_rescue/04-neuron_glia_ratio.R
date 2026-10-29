# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

source(here("code/analysis/18_rescue/rescue_plot_utils.R"))


# Config ----
fld_out <- here("processed-data/rds/18_rescue")

# Same threshold as 03-donor_spd_RESCUE_composition.R: donor x domain
# combinations with fewer spots are excluded from the across-donor summary
min_spots <- 10

# Broad classes making up each compartment (EndoMural and Unattributed are in
# neither)
neuron_classes <- c("Excit", "Inhib")
glia_classes <- c("Astro", "Oligo", "OPC", "Micro")


# Load data ----
## Slim key-indexed RESCUE proportions, from 01-merge_RESCUE_to_spe.R ----
prop_df <- readRDS(
  here("processed-data/rds/18_rescue", "RESCUE_prop_df.rds")
)

ct_cols <- setdiff(
  grep("^RESCUE_", colnames(prop_df), value = TRUE),
  c("RESCUE_nUMI", "RESCUE_unattributed_ratio")
)

# error prevention
stopifnot(!anyNA(prop_df[, ct_cols]))
stopifnot(all(c(neuron_classes, glia_classes) %in% rescue_broad_class(ct_cols)))


# Spot-level neuronal and glial proportions ----------------------------------
neuron_cols <- ct_cols[rescue_broad_class(ct_cols) %in% neuron_classes]
glia_cols <- ct_cols[rescue_broad_class(ct_cols) %in% glia_classes]

cat("Neuronal cell types:", str_remove(neuron_cols, "^RESCUE_"), "\n")
cat("Glial cell types:", str_remove(glia_cols, "^RESCUE_"), "\n")

spot_df <- prop_df |>
  mutate(
    neuron_prop = rowSums(across(all_of(neuron_cols))),
    glia_prop = rowSums(across(all_of(glia_cols)))
  ) |>
  select(key, brnum, dx, spd_label, neuron_prop, glia_prop)


# Donor x spatial domain neuron-to-glia ratio --------------------------------
## Ratio of mean proportions within each donor x domain (and donor overall):
##   ratio = mean(neuron_prop) / mean(glia_prop)
## A ratio of means rather than a mean of per-spot ratios, since many spots
## have near-zero glial proportion (per-spot ratios are unstable / infinite).
## Unattributed expression cancels in the ratio, so it is identical with or
## without renormalizing to attributed signal.
summarize_ratio <- function(df) {
  df |>
    summarize(
      n_spots = n(),
      mean_neuron_prop = mean(neuron_prop),
      mean_glia_prop = mean(glia_prop),
      .groups = "drop"
    ) |>
    mutate(neuron_glia_ratio = mean_neuron_prop / mean_glia_prop)
}

donor_spd_ratio <- bind_rows(
  spot_df |>
    group_by(brnum, dx, spd_label) |>
    summarize_ratio() |>
    mutate(spd_label = as.character(spd_label)),
  spot_df |>
    group_by(brnum, dx) |>
    summarize_ratio() |>
    mutate(spd_label = "All spots")
) |>
  mutate(
    spd_label = factor(spd_label, levels = c(levels(prop_df$spd_label), "All spots")),
    dx = factor(dx, levels = c("ntc", "scz")),
    pass_min_spots = n_spots >= min_spots
  ) |>
  arrange(spd_label, dx, brnum)

# error prevention
stopifnot(all(donor_spd_ratio$mean_glia_prop > 0))
stopifnot(!anyDuplicated(donor_spd_ratio[, c("brnum", "spd_label")]))

cat(sprintf(
  "Donor x domain combinations with < %d spots (excluded from summary): %d\n",
  min_spots, sum(!donor_spd_ratio$pass_min_spots)
))
print(donor_spd_ratio |> filter(!pass_min_spots), n = Inf)


# Median and IQR across donors per spatial domain -----------------------------
spd_ratio_summary <- donor_spd_ratio |>
  filter(pass_min_spots) |>
  group_by(spd_label) |>
  summarize(
    n_donors = n(),
    median_ratio = median(neuron_glia_ratio),
    q1_ratio = quantile(neuron_glia_ratio, 0.25),
    q3_ratio = quantile(neuron_glia_ratio, 0.75),
    iqr_ratio = IQR(neuron_glia_ratio),
    .groups = "drop"
  )

print(spd_ratio_summary, width = Inf)


# Save -----------------------------------------------------------------------
write_csv(
  donor_spd_ratio,
  file.path(fld_out, "RESCUE_neuron_glia_ratio_donor_spd.csv")
)

write_csv(
  spd_ratio_summary,
  file.path(fld_out, "RESCUE_neuron_glia_ratio_spd_summary.csv")
)

print("Finished computing RESCUE neuron-to-glia ratios per donor and spatial domain")


# Session info ----
sessioninfo::session_info()
