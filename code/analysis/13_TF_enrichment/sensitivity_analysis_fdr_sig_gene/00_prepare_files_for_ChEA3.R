# Sensitivity analysis: same as ../00_prepare_files_for_ChEA3.R, but using
# layer-restricted DEGs with FDR < 0.10 instead of nominal p-value < 0.05

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Minimum number of DEGs needed to run ChEA3 for a domain-direction
# NOTE: many domains have no or very few DEGs at FDR < 0.10. TF ranks from
# 1-2 genes are not interpretable, so these gene sets are skipped.
min_n_genes <- 5

out_folder <- here(
  "processed-data/rds/13_TF_enrichment/sensitivity_analysis_fdr_sig_gene",
  "input"
)
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

# Save files ----
## Load layer-restricted DE results ----
spd_deg_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
)

fdr_10_df <- spd_deg_df |>
  filter(adj.P.Val < 0.10) |>
  mutate(direction = if_else(logFC > 0, "Up", "Down"))

# Number of DEGs per domain and direction
fdr_10_df |>
  count(PRECAST_spd, direction) |>
  complete(
    PRECAST_spd = unique(spd_deg_df$PRECAST_spd),
    direction = c("Up", "Down"),
    fill = list(n = 0)
  ) |>
  mutate(run_ChEA3 = n >= min_n_genes) |>
  print(n = Inf)

fdr_10_df |>
  group_by(PRECAST_spd, direction) |>
  filter(n() >= min_n_genes) |>
  group_walk(
    .f = function(.x, .y) {
      .x |>
        pull(gene) |>
        write_lines(
          file.path(
            out_folder,
            paste0(
              "ChEA3_input-",
              # Keep the same naming as the main analysis, e.g. SpD07_L1
              gsub("[^A-Za-z0-9_]", "_", str_remove(.y$PRECAST_spd, "/M$")),
              "-", .y$direction,
              ".txt"
            )
          )
        )
    }
  )

# Session info ----
session_info()
