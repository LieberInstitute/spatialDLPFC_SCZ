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
        "SpD07-L1/M",
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
  pattern = "layer_restricted_logFC_.*\\.csv",
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
          "(?<=layer_restricted_logFC_).*?(?=\\.csv)"
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
    "layer_specific_genes_nom_p.csv"
  )
) |> mutate(layer_specific_nom_p_005 = TRUE)

unique_genes_fdr <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_specific_genes_fdr_010.csv"
  )
) |> mutate(layer_specific_fdr_010 = TRUE)

spd_deg_df <- spd_deg_df |>
  left_join(
    unique_genes_nom,
    by = c("gene_id" = "ENSEMBL", "PRECAST_spd" = "spd")
  ) |>
  left_join(
    unique_genes_fdr,
    by = c("gene_id" = "ENSEMBL", "PRECAST_spd" = "spd")
  )

#|>
spd_deg_df[
  which(is.na(spd_deg_df$layer_specific_nom_p_005)),
  "layer_specific_nom_p_005"
] <- FALSE
spd_deg_df[
  which(is.na(spd_deg_df$layer_specific_fdr_010)),
  "layer_specific_fdr_010"
] <- FALSE
# mutate(
#   layer_specific_nom_p_005 = if_else(is.na(layer_specific), FALSE, TRUE)
# )

sum(spd_deg_df$layer_specific_nom_p_005)
# Suppose to be 4073
sum(spd_deg_df$layer_specific_fdr_010)
# Suppose to be 913

sum(!(spd_deg_df$layer_specific_nom_p_005) &
  spd_deg_df$layer_specific_fdr_010)


# TO NOTE:
# It's possibel that the layer-specific DEG by FDR is not a layer-specific DEG by nom p.
# THis is because there could be a gene that are sig at nom-p for all the layers, but fdr sig only in one layer.
# Taking MAPK3 as an example, see code below:

# spd_deg_df |> filter(
#   spd_deg_df$layer_specific_fdr_010,
#   !(spd_deg_df$layer_specific_nom_p_005)
# ) |> View()
# spd_deg_df |> filter(gene == "MAPK3") |> View()



## Save to csv ----
write_csv(
  spd_deg_df,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  ),
  na = "",
  quote = "all"
)

# Session info ----
session_info()
