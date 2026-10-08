# Run ChEA3 via its API (instead of manually through the web portal)
# NOTE: the API returns the same Integrated--meanRank table as the portal's
#  "Integrated_meanRank" download, e.g. identical ranks for SpD01_WMtz-Up in
#  the main analysis.
#  See https://maayanlab.cloud/chea3/index.html#content4-z

# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(sessioninfo)
})

chea_folder <- here(
  "processed-data/rds/13_TF_enrichment/sensitivity_analysis_fdr_sig_gene"
)

# Run ChEA3 ----
input_df <- data.frame(
  full_path = list.files(
    file.path(chea_folder, "input"),
    pattern = "ChEA3_input-.*\\.txt",
    full.names = TRUE
  )
) |>
  mutate(
    spd = str_split_i(basename(full_path), "-|\\.txt", 2),
    direction = str_split_i(basename(full_path), "-|\\.txt", 3) |>
      str_to_lower()
  )

input_df |>
  pwalk(
    .f = function(full_path, spd, direction) {
      genes <- read_lines(full_path)

      response <- POST(
        url = "https://maayanlab.cloud/chea3/api/enrich/",
        body = list(query_name = "gene_set_query", gene_set = genes),
        encode = "json"
      )
      stop_for_status(response)

      results <- content(response, "text", encoding = "UTF-8") |>
        fromJSON()

      results[["Integrated--meanRank"]] |>
        write_tsv(
          file.path(
            chea_folder,
            paste0("Integrated_meanRank-", spd, "-", direction, ".tsv")
          )
        )
    }
  )

# Session info ----
session_info()
