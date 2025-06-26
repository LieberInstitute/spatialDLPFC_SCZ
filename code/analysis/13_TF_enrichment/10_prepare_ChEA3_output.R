# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load data ----
chea_folder <- here("processed-data/rds/13_TF_enrichment")

file_df <- data.frame(
  full_path = list.files(
    chea_folder,
    pattern = "Integrated_meanRank-.*\\.tsv",
    full.names = TRUE
  )
) |>
  mutate(
    rel_path = str_remove(full_path, paste0(chea_folder, "/")),
    spd = str_split_i(rel_path, "-|\\.tsv", 2),
    direction = str_split_i(rel_path, "-|\\.tsv", 3)
    |> str_to_lower()
  )

raw_chea_df <- file_df |>
  pmap(
    .f = function(full_path, rel_path, spd, direction) {
      read_tsv(full_path, show_col_types = FALSE) |>
        mutate(
          spd = spd,
          direction = direction
        ) |>
        select(spd, direction, everything()) |>
        select(-`Query Name`)
    }
  ) |>
  bind_rows()

# Write output ----
raw_chea_df |>
  write_csv(
    here("processed-data/rds/13_TF_enrichment", "ChEA3_result_compiled.csv")
  )


# Load session info ----
session_info()
