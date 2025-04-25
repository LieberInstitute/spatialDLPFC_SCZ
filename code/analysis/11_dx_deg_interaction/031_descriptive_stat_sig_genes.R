# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggrepel)
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

##


# Descriptive Statistics ----
## Nomial p-value < 0.05 ----
nom_p_vec <- spd_deg_list |>
  map_int(
    ~ sum(.x$P.Value < 0.05, na.rm = TRUE)
  )

# SpD01-WMtz SpD02-L3/4   SpD03-L6   SpD04-WM   SpD05-L5 SpD06-L2/3
#       1581        284        283       2521        468        782
#   SpD07-L1
#       1059


## FDR < 0.1 ----
fdr_10_vec <- spd_deg_list |>
  map_int(
    ~ sum(.x$adj.P.Val < 0.1, na.rm = TRUE)
  )

# SpD01-WMtz SpD02-L3/4   SpD03-L6   SpD04-WM   SpD05-L5 SpD06-L2/3
#         30          0          0        894          1          2
#   SpD07-L1
#         17


cbind(
  nom_p_vec,
  fdr_10_vec
)

#            nom_p_vec fdr_10_vec
# SpD01-WMtz      1581         30
# SpD02-L3/4       284          0
# SpD03-L6         283          0
# SpD04-WM        2521        894
# SpD05-L5         468          1
# SpD06-L2/3       782          2
# SpD07-L1        1059         17


# Upset plots ----
## Nominal p-value < 0.05 ----
nom_p_geneList <- spd_deg_list |>
  map(
    ~ .x$gene_id[.x$P.Value < 0.05]
  )

pdf(
  here(
    "plots/11_dx_deg_interaction",
    "upset_layer_gene_nominal_p_value.pdf"
  ),
  width = 10,
  height = 8
)
upset(
  fromList(nom_p_geneList),
  nsets = 7,
  # sets = sort(
  #   names(enrich_list),
  #   decreasing = TRUE
  # ),
  keep.order = TRUE,
  nintersects = NA,
  order.by = c("freq"),
  # main.bar.color = "blue",
  # sets.bar.color = "red",
  # mainbar.y.label = "Enriched Genes",
  sets.x.label = "Spatial Domains"
)
dev.off()

## FDR p-value < 0.10 ----

fdr_10_geneList <- spd_deg_list |>
  map(
    ~ .x$gene_id[.x$adj.P.Val < 0.10]
  )


pdf(
  here(
    "plots/11_dx_deg_interaction",
    "upset_layer_gene_fdr_10.pdf"
  ),
  width = 10,
  height = 8
)
upset(
  fromList(fdr_10_geneList),
  nsets = 7,
  # sets = sort(
  #   names(enrich_list),
  #   decreasing = TRUE
  # ),
  keep.order = TRUE,
  nintersects = NA,
  order.by = c("freq"),
  # main.bar.color = "blue",
  # sets.bar.color = "red",
  # mainbar.y.label = "Enriched Genes",
  sets.x.label = "Spatial Domains"
)
dev.off()


# Idenitfy Layer specific genes ----
# Nom p-value ----
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

write_csv(
  unique_genes_nom,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_uniquely_specific_genes_nom_p.csv"
  )
)

# FDR p-value 10 ----

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
