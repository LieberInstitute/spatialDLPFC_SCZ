# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(CellChat)
  library(here)
  library(sessioninfo)
})

## Source in-house function ----
source(here(
  "code/analysis/Layer_layer_communication/",
  "vulnerablility_charact/fun_extract_pathway_diff_mat.r"
))

# Load CellChat objects ----
ntc_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "ntc_cellchat_annotated.rds"
  )
)

scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_annotated.rds"
  )
)

# Find all pathways ----
cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "merged_cellchat_annotated.rds"
  )
)

cellchat_pathways <- rankNet(cellchat, return.data = TRUE)$signaling.contribution

## Characterize pathways ----
nrow(cellchat_pathways)
# [1] 60
cellchat_diff_pathways <- unique(cellchat_pathways$name) |> as.character()

# Extract each pathways ----
df_all_pathways <- cellchat_diff_pathways |>
  # iterate over all pathways
  map_dfr(.f = function(.path) {
    # Pull the pathway diff matrix
    ret_raw <- extract_diff_pathway(
      ntc_cellchat, scz_cellchat,
      signaling = .path, measure = "weight",
      thresh = 0.05
    )

    # format the pathway
    ret_df_long <- ret_raw |>
      as.data.frame(check.names = FALSE) |>
      rownames_to_column(var = "sender") |>
      pivot_longer(
        cols = -sender,
        names_to = "receiver",
        values_to = "diff_value"
      ) |>
      # Add pathway information
      mutate(
        pathway = .path
      ) |>
      # Reorder columns to have pathway first
      select(pathway, sender, receiver, diff_value)

    return(ret_df_long)
  })

# error prevention
stopifnot(
  nrow(df_all_pathways) ==
    length(cellchat_diff_pathways) * 7 * 7 # 7 SpDs
)

# Save csv_df ----
df_all_pathways |>
  write_csv(
    here(
      "code/analysis/Layer_layer_communication/vulnerablility_charact",
      "diff_value_pathway_all_pairs.csv"
    ),
    quote = "all"
  )

# Session info ----
session_info()
