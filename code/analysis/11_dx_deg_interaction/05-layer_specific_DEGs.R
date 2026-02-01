# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

# Load data ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_restricted_logFC_.*\\.csv",
  full.names = TRUE
)

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

spd_deg_list <-
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
  map(
    ~ read_csv(.x)
  )

# Idenitfy Layer specific genes ----
## Nom p-value ----

nom_p_geneList <- spd_deg_list |>
  map(
    ~ .x$gene_id[.x$P.Value < 0.05]
  )

unique_genes_nom <- map(
  names(nom_p_geneList),
  ~ setdiff(
    nom_p_geneList[[.x]],
    # Reduce(union, nom_p_geneList[names(nom_p_geneList) != .x])
    unlist(nom_p_geneList[names(nom_p_geneList) != .x]) |> unique()
  )
) |>
  set_names(names(nom_p_geneList))


unique_genes_nom <- imap_dfr(
  unique_genes_nom,
  .f = function(.x, idx) {
    # browser()
    # Convert genes from ensembl to symbol
    # NOTE: some ensembl ID may not be converted.
    # clusterProfiler::bitr(
    #   .x,
    #   fromType = "ENSEMBL", toType = "SYMBOL",
    #   OrgDb = "org.Hs.eg.db",
    #   drop = FALSE
    # ) |>
    data.frame(
      ENSEMBL = .x,
      spd = idx
    )
  }
)

length(unique_genes_nom$ENSEMBL)
# [1] 4073

stopifnot(length(unique(unique_genes_nom$ENSEMBL)) == nrow(unique_genes_nom))
stopifnot(all(!duplicated(unique_genes_nom$ENSEMBL)))

write_csv(
  unique_genes_nom,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_specific_genes_nom_p.csv"
  ),
  quote = "all"
)

## FDR p-value 10 ----
fdr_p_geneList <- spd_deg_list |>
  map(
    ~ .x$gene_id[.x$adj.P.Val < 0.10]
  )

unique_genes_fdr <- map(
  names(fdr_p_geneList),
  ~ setdiff(
    fdr_p_geneList[[.x]],
    # Reduce(union, nom_p_geneList[names(nom_p_geneList) != .x])
    unlist(fdr_p_geneList[names(fdr_p_geneList) != .x]) |> unique()
  )
) |>
  set_names(names(fdr_p_geneList))

unique_genes_fdr <- imap_dfr(
  unique_genes_fdr,
  .f = function(.x, idx) {
    # browser()
    # Convert genes from ensembl to symbol
    # NOTE: some ensembl ID may not be converted.
    # clusterProfiler::bitr(
    #   .x,
    #   fromType = "ENSEMBL", toType = "SYMBOL",
    #   OrgDb = "org.Hs.eg.db",
    #   drop = FALSE
    # ) |>
    if (length(.x) != 0) {
      data.frame(
        ENSEMBL = .x,
        spd = idx
      )
    } else {
      data.frame(
        ENSEMBL = character(),
        spd = character()
      )
    }
  }
)

length(unique_genes_fdr$ENSEMBL)
# [1] 913

stopifnot(length(unique(unique_genes_fdr$ENSEMBL)) == nrow(unique_genes_fdr))
stopifnot(all(!duplicated(unique_genes_fdr$ENSEMBL)))

write_csv(
  unique_genes_fdr,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_specific_genes_fdr_010.csv"
  ),
  quote = "all"
)


# unique_genes_fdr <- map(
#   names(fdr_10_geneList),
#   ~ setdiff(
#     fdr_10_geneList[[.x]],
#     # Reduce(union, fdr_10_geneList[names(fdr_10_geneList) != .x])
#     unlist(fdr_10_geneList[names(fdr_10_geneList) != .x]) |> unique()
#   )
# ) |>
#   set_names(names(fdr_10_geneList))


# unique_genes_fdr <- imap_dfr(
#   unique_genes_fdr,
#   .f = function(.x, idx) {
#     # browser()
#     # Convert genes from ensembl to symbol
#     # NOTE: some ensembl ID may not be converted.
#     data.frame(
#       ENSEMBL = .x,
#       spd = idx
#     )
#   }
# )

# write_csv(
#   unique_genes_fdr,
#   here(
#     "processed-data/rds/11_dx_deg_interaction",
#     "layer_specific_genes_fdr.csv"
#   )
# )

# Session Info ----
session_info()
