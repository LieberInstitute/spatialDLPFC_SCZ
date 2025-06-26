# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Save files ----
## Load layer-restricted DE results ----
spd_deg_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
)

unique(spd_deg_df$PRECAST_spd) |>
  walk(
    .f = function(.spd) {
      nom_05_df <- spd_deg_df |>
        filter(PRECAST_spd == .spd) |>
        filter(P.Value < 0.05)

      nom_05_df |>
        filter(logFC > 0) |>
        pull(gene) |>
        write_lines(
          here(
            "processed-data/rds/13_TF_enrichment/input",
            paste0(
              "ChEA3_input-",
              gsub("[^A-Za-z0-9_]", "_", .spd),
              "-Up.txt"
            )
          )
        )

      nom_05_df |>
        filter(logFC < 0) |>
        pull(gene) |>
        write_lines(
          here(
            "processed-data/rds/13_TF_enrichment/input",
            paste0(
              "ChEA3_input-",
              gsub("[^A-Za-z0-9_]", "_", .spd),
              "-Down.txt"
            )
          )
        )
    }
  )


# Create .txt file ----



# Session info ----
session_info()
