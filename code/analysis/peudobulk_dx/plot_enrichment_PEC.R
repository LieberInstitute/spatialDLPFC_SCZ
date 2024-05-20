# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(pheatmap)
  library(sessioninfo)
})

# Read Data ----
## Dx DEG genes ----
gene_df <- read_csv(
  #   here(
  #     "processed-data/PB_dx_genes/",
  #     "test_PRECAST_07.csv"
  #   )
  # )

  "~/Downloads/test_PRECAST_07.csv"
)

sig_gene_df <- gene_df |>
  filter(fdr_ntc <= 0.05)

## PEC SpD Enrichment ----
# TODO: find the raw dataset for all genes instead of some genes
# pec_spd_df <- read.csv(
# here(
# "code/analysis/visium_spatial_clustering",
# "TableS8_sig_genes_FDR5perc_enrichment.csv"

# )
# Local path
load(
  "~/modeling_results_BayesSpace_k09.Rdata"
)

pec_spd_df <- modeling_results$enrichment # |> str()
# filter(spatial_domain_resolution == "Sp09") #


# pec_spd_gene_sub <- pec_spd_df |>
#   filter(pec_spd_df$ensembl %in% sig_gene_df$ensembl)


# Heatmap
# Orientation: gene (row) by spd (columns) with fill (-log10(p), maybe a directionality of logFC)

heatmap_pec_spd_df <- pec_spd_df |>
  mutate(
    across(
      starts_with("p_value"),
      ~ -1 * log10(.x),
      .names = "neg_log10_{.col}"
    )
  ) |>
  filter(ensembl %in% sig_gene_df$ensembl) |>
  select(starts_with("neg_log10_"))

# stopifnot()

# mutate(
#   neg_log10_p = -1 * log10(pval)
# ) |>
# pivot_wider(
#   id_cols = c("ensembl"),
#   names_from = "test",
#   values_from = "neg_log10_p"
# ) |>
# column_to_rownames(var = "ensembl") |>
# select(starts_with("Sp09")) |>

heatmap_pec_spd_df |>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cellwidth = 5,
    cellheight = 5
  )




## PEC snRNA Enrichment ----

# /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
# here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_velm_registration.Rdata")


load(
  # TODO: this is the wrong dataset, needs update after confirming with Louise
  here(
    "sn_velm_registration.Rdata"
  )
)

pec_snRNA_df <- sn_hc_registration$enrichment

pec_snRNA_df |>
  mutate(
    across(
      starts_with("p_value"),
      ~ -1 * log10(.x),
      .names = "neg_log10_{.col}"
    )
  ) |>
  filter(ensembl %in% sig_gene_df$ensembl) |>
  select(starts_with("neg_log10_"))|>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cellwidth = 5,
    cellheight = 5
  )


# Session Info ----
session_info()