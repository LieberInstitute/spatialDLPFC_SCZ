# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(spatialLIBD)
  library(sessioninfo)
})

# Load Data ----
## load Layer_specific DEGs ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

names(spd_files) <- str_extract(
  spd_files,
  "(?<=layer_specific_logFC_).*?(?=\\.csv)"
)

spd_deg_list <- map(
  spd_files,
  ~ read_csv(.x, col_types = cols()) |>
    mutate(gene = AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )) |>
    filter(
      P.Value < 0.05
    )
)

spd_deg_list |> map(~ .x |> nrow())

## Load Huuki-Myers cell type -markers ----
load(
  here(
    "code/analysis/12_deg_integration",
    "huuki_myers_cell_type_enrichment_cellType_hc.Rdata"
  )
)

cell_type_enrich_df <- results_specificity

# Enrichment Test ----
res <- spatialLIBD::gene_set_enrichment(
  gene_list = spd_deg_list |> map(~ .x |> pull(gene_id)),
  modeling_results = list(
    "enrichment" = cell_type_enrich_df
  ),
  model_type = "enrichment",
  fdr_cut = 0.05
)

gene_set_enrichment_plot(
  res,
  PThresh = 12,
  ORcut = 3,
  enrichOnly = FALSE,
  cex = 1.5 # control the size of the text
)
