# Load Library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  library(here)
})

# Load Data ----
## Load DEGs ----

load(
  here(
    "processed-data/BrainSeqV2",
    "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda"
  )
)
# NOTE: contains outGene (what we need w. qSV) and outGene0 and outGeneNoAdjust
# https://github.com/LieberInstitute/brainseq_phase2?tab=readme-ov-file#dxstats_dlpfc_filtered_qsva_genelevel_nohgoldqsv_matchdlpfcrda

gene_df <-  outGene |> filter(
  adj.P.Val < 0.10
)




## Load SPG enrichment res ----
### Neuropil ----
neuropil_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "neuropil_dx_res.csv"
  )
)

### Neun ----
neun_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "neun_dx_res.csv"
  )
)

### Vasculature ----
vasc_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "vasc_dx_res.csv"
  )
)

### pnn ----
pnn_df <- read_csv(
  here(
    "processed-data/image_processing/enrichment",
    "pnn_dx_res.csv"
  )
)

# Exact test for enrichment ----
## Neuropil ----
neuropil_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `BrainSeqV2` = gene_df  |> pull(ensemblID),
    `upreg genes` = gene_df |> filter(logFC > 0) |> pull(ensemblID),
    `downreg genes` = gene_df |> filter(logFC < 0) |> pull(ensemblID)
  ),
  modeling_results = list("enrichment" = neuropil_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "Neuropil"
  ) #|>
  # gene_set_enrichment_plot(
  #   PThresh = 12,
  #   ORcut = 3,
  #   enrichOnly = FALSE,
  #   cex = 1.5 # control the size of the text
  # ) +
  # title("SCZ-DEGs enriched in Neuropil")

## Neun ----
neun_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `BrainSeqV2` = gene_df  |> pull(ensemblID),
    `upreg genes` = gene_df |> filter(logFC > 0) |> pull(ensemblID),
    `downreg genes` = gene_df |> filter(logFC < 0) |> pull(ensemblID)
  ),
  modeling_results = list("enrichment" = neun_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "Neun"
  )

## Vasculature ----
vasc_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `BrainSeqV2` = gene_df  |> pull(ensemblID),
    `upreg genes` = gene_df |> filter(logFC > 0) |> pull(ensemblID),
    `downreg genes` = gene_df |> filter(logFC < 0) |> pull(ensemblID)
  ),
  modeling_results = list("enrichment" = vasc_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "Vasculature"
  )


## PNN ----
pnn_enrich <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `BrainSeqV2` = gene_df  |> pull(ensemblID),
    `upreg genes` = gene_df |> filter(logFC > 0) |> pull(ensemblID),
    `downreg genes` = gene_df |> filter(logFC < 0) |> pull(ensemblID)
  ),
  modeling_results = list("enrichment" = pnn_df),
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  filter(
    test == TRUE
  ) |>
  mutate(
    test = "PNN"
  )

# Combine results ----
bind_rows(
  neuropil_enrich,
  neun_enrich,
  vasc_enrich,
  pnn_enrich
) |>
  # mutate(
  #   test = factor(
  #     test,
  #     levels = c("Neuropil", "Neun", "Vasculature", "PNN")
  #   )
  # ) |>
  gene_set_enrichment_plot(
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    cex = 1.5 # control the size of the text
  ) +
  title("BrainSeqV2 dx-DEGs enriched in SPG-called categories")