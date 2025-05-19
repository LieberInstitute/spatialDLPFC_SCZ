# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load data ----
## Load spd annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(
    anno_lab = factor(
      paste0(gsub("spd", "SpD", spd), "-", label),
      levels = c(
        "SpD07-L1",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      )
    )
  )


## Load gene_names ----
# NOTE: to connect ensemble ID to gene symbol
adj_deg_df_raw <- read_csv(here(
  "processed-data/rds/10_dx_deg_adjust_spd",
  "dx-deg_PRECAST07.csv"
))

## Load layer-specific DEGs files ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

spd_deg_df <-
  spd_files |>
  set_names(
    # Retrieve annotated names
    nm = spd_anno_df[
      match(
        str_extract(
          spd_files,
          "(?<=layer_specific_logFC_).*?(?=\\.csv)"
        ),
        spd_anno_df$spd
      ),
      "anno_lab",
      drop = TRUE
    ]
  ) |>
  imap_dfr(
    ~ read_csv(.x) |>
      right_join(
        adj_deg_df_raw |> select(ensembl, gene),
        by = c("gene_id" = "ensembl")
      ) |>
      mutate(PRECAST_spd = .y)
  )

# Error prevention: genes doesn't have gene symbols
stopifnot(spd_deg_df |> filter(is.na(gene)) == 0)

# Annotate Layer-specific DEGs ----
unique_genes_nom <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_uniquely_specific_genes_nom_p.csv"
  )
)

spd_deg_df <- spd_deg_df |>
  left_join(
    unique_genes_nom |> select(-SYMBOL) |> mutate(layer_specific = TRUE),
    by = c("gene_id" = "ENSEMBL", "PRECAST_spd" = "spd")
  ) |>
  mutate(layer_specific = if_else(is.na(layer_specific), FALSE, TRUE))

sum(spd_deg_df$layer_specific)

## Save to csv ----
write_csv(
  spd_deg_df,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  ),
  na = ""
)

# Session info ----
session_info()
