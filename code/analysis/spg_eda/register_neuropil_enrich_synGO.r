# NOTE: imput data is Maddy's enrichment file

# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
  library(sessioninfo)
})

# Load Data ----
## Load Maddy's enrichment res ----
neuropil_enrich_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "neuropil_dx_res.csv"
  )
)


## Load SynGO Genes ----
synGO_genes <- read_excel(
  here::here(
    "processed-data/SynGO/syngo_genes.xlsx"
  )
)

# Finding overlapping genes ---
## Overlapping in Sig Genes ----
sig_synGO_gene_df <- neuropil_enrich_df |>
  filter(
    fdr_TRUE < 0.05
  ) |>
  inner_join(
    synGO_genes,
    by = c("ensembl" = "ensembl_id")
  )


sig_synGO_gene_df$logFC_TRUE |> hist()
sig_synGO_gene_df$t_stat_TRUE |> hist()


## Overlapping in backgrounds ----
# This looks fine. like not drastically bad.
hist(neuropil_enrich_df$p_value_TRUE, breaks = 100)

## Plots ----
# Create Venn Diagram ----
library(VennDiagram)

# Define the sets
sig_genes <- neuropil_enrich_df |>
  filter(fdr_TRUE < 0.05) |>
  pull(ensembl)

background_genes <- neuropil_enrich_df |>
  pull(ensembl)

syngo_genes <- synGO_genes |>
  pull(ensembl_id)

# Generate Venn Diagram
venn.diagram(
  x = list(
    "Significant Genes" = sig_genes,
    "All Genes in PB data" = background_genes,
    "SynGO Genes" = syngo_genes
  ),
  filename = NULL,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Venn Diagram of Gene Overlap",
  disable.logging = TRUE
)


# Session Info ----
sessioninfo::session_info()
