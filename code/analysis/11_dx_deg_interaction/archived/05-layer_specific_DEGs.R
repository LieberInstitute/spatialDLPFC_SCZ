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
  pattern = "layer_specific_logFC_.*\\.csv",
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
          "(?<=layer_specific_logFC_).*?(?=\\.csv)"
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
    Reduce(union, nom_p_geneList[names(nom_p_geneList) != .x])
  )
) |>
  set_names(names(nom_p_geneList))


unique_genes_nom <- imap_dfr(
  unique_genes_nom,
  .f = function(.x, idx) {
    # browser()
    # Convert genes from ensembl to symbol
    # NOTE: some ensembl ID may not be converted.
    clusterProfiler::bitr(
      .x,
      fromType = "ENSEMBL", toType = "SYMBOL",
      OrgDb = "org.Hs.eg.db",
      drop = FALSE
    ) |>
      mutate(
        spd = idx
      )
  }
)

length(unique_genes_nom$ENSEMBL)
# [1] 4114

log_df <- spd_deg_list |>
  imap_dfr(
    ~ .x |>
      filter(P.Value < 0.05) |>
      filter(gene_id %in% unique_genes_nom$ENSEMBL) |>
      mutate(spd = .y)
  )
stopifnot(nrow(log_df) == length(unique_genes_nom$ENSEMBL))

write_csv(
  unique_genes_nom,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_uniquely_specific_genes_nom_p.csv"
  )
)

## FDR p-value 10 ----

unique_genes_fdr <- map(
  names(fdr_10_geneList),
  ~ setdiff(
    fdr_10_geneList[[.x]],
    Reduce(union, fdr_10_geneList[names(fdr_10_geneList) != .x])
  )
) |>
  set_names(names(fdr_10_geneList))


unique_genes_fdr <- imap_dfr(
  unique_genes_fdr,
  .f = function(.x, idx) {
    # browser()
    # Convert genes from ensembl to symbol
    # NOTE: some ensembl ID may not be converted.
    clusterProfiler::bitr(
      .x,
      fromType = "ENSEMBL", toType = "SYMBOL",
      OrgDb = "org.Hs.eg.db",
      drop = FALSE
    ) |>
      mutate(
        spd = idx
      )
  }
)

write_csv(
  unique_genes_fdr,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_uniquely_specific_genes_fdr.csv"
  )
)

# Session Info ----
session_info()
