# Shared helpers for RESCUE visualization scripts (sourced, not run directly)

# Broad cell class of a RESCUE cell type column, e.g.
# "RESCUE_Excit_03" -> "Excit", "RESCUE_Unattributed" -> "Unattributed"
rescue_broad_class <- function(ct_col) {
  str_remove(ct_col, "^RESCUE_") |> str_remove("_\\d+$")
}

# Stacking order of broad classes (neurons first, unattributed last)
broad_class_levels <- c(
  "Excit", "Inhib", "Astro", "Oligo", "OPC", "Micro", "EndoMural",
  "Unattributed"
)

# Broad-class palette for compositional plots; "Unattributed" is a neutral
# gray since it is not a cell type.
# NOTE: Inhib and Oligo are similar hues; `broad_class_levels` keeps them
# non-adjacent in stacked bars.
broad_class_colors <- c(
  Astro = "#4DAF6E", # green
  EndoMural = "#E7609E", # pink
  Micro = "#5E3C99", # purple
  Oligo = "#E07B39", # orange
  OPC = "#D4B243", # mustard
  Excit = "#2E6DB4", # blue (all excitatory types grouped)
  Inhib = "#E05A42", # red-orange
  Unattributed = "#b5b4ae" # neutral gray
)

# Conditional composition excluding unattributed expression: drop the
# "Unattributed" class and renormalize the remaining classes to sum to 1
# within each group. Applied to already-averaged proportions (`mean_prop`),
# i.e. the share of *attributed* signal from each cell type, which equals
#   mean(p_ct) / mean(1 - p_unattributed)
# over the spots in the group (avoids inflating spots with almost no
# attributed signal, as per-spot renormalization would).
to_conditional_comp <- function(df, group_vars) {
  df |>
    filter(broad_class != "Unattributed") |>
    group_by(across(all_of(group_vars))) |>
    mutate(mean_prop = mean_prop / sum(mean_prop)) |>
    ungroup() |>
    mutate(broad_class = fct_drop(broad_class))
}

# Shared fill scale for broad-class stacked bars
scale_fill_broad_class <- function() {
  scale_fill_manual(
    values = broad_class_colors, name = "Broad class",
    breaks = broad_class_levels
  )
}

# Two versions of every compositional plot: as reported by RESCUE (incl.
# unattributed) and conditional on attributed signal (excl. unattributed)
comp_versions <- tribble(
  ~version, ~file_suffix, ~y_label,
  "with_unattributed", "", "Mean RESCUE proportion",
  "excl_unattributed", "_excl_unattributed", "Mean conditional proportion\n(excl. unattributed)"
)

dx_colors <- c(ntc = "blue", scz = "red")
dx_labels <- c(ntc = "NTC", scz = "SCZ")

# Representative samples: Br8667 (NTC) and Br5973 (SCZ)
# (same pair used throughout the paper, e.g. 03_visium_spatial_clustering)
rep_sample_id <- c("V13M06-342_D1", "V13M06-343_D1")
