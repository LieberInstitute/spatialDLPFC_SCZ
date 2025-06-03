# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(readxl)
  library(sessioninfo)
})

# Load Data ----
## Load Huuki-Myer et al. (2024) cell-type markers ----
cell_type_enrich_df <- readRDS(
  # NOTE: if run on JHPCE
  # file.path(
  #   # PEC Study folder
  #   "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/",
  #   "processed-data/rdata/spe/14_spatial_registration_PEC",
  #   "registration_stats_LIBD.rds"
  # )
  # NOTE: if run on Boyi's local machine
  here(
    "code/analysis/12_cross-study_enrichment",
    "registration_stats_LIBD.rds"
  )
)

## Xenium panel data ----
genes_xenium <- read_xlsx(
  here(
    "code/xenium_panel_design", "Xenium_SHK_celltype_REannot_2025-04-13.xlsx"
  ),
  sheet = "combined_Xenium_SCZ_ProbeSelect"
) |> select(
  -`...1`
)



# Create DF for L6 related marker genes ----
## Retrieve # of marker genes per cell type ----
tstats <- cell_type_enrich_df[, grep("[f|t]_stat_", colnames(cell_type_enrich_df))]
colnames(tstats) <- gsub("[f|t]_stat_", "", colnames(tstats))
fdrs <- cell_type_enrich_df[, grep("fdr_", colnames(cell_type_enrich_df))]
colnames(fdrs) <- gsub("fdr_", "", colnames(fdrs))

## Identify Cell types with significant markers ----
all_cell_types_df <- colnames(tstats) |>
  set_names() |>
  map(
    ~ cell_type_enrich_df[
      # Rows
      which(tstats[, .x] > 0 & fdrs[, .x] < 0.05),
      # Columns
      c("ensembl", "gene")
    ] |> mutate(
      PEC_cell_type = .x
    )
  ) |>
  bind_rows() |>
  group_by(gene) |>
  summarize(
    ensembl = first(ensembl),
    n_PEC_cell_type = n(),
    PEC_cell_type = paste(unique(PEC_cell_type), collapse = ", "),
    .groups = "drop"
  )

left_join(
  genes_xenium,
  all_cell_types_df,
  by = c("Ensembl ID" = "ensembl")
) |>
  # filter(!is.na(n_PEC_cell_type)) |>
  select(
    `Ensembl ID`, Gene, n_PEC_cell_type, PEC_cell_type, cell_type_updated, layer_marker, Note, everything()
  ) |>
  select(
    - `gene`,
    "Anno_layer_marker" = "layer_marker",
    "Anno_cell_type" = "cell_type_updated",
    "Anno_cell_type_obselete" = "cell_type"
  ) |> 
write_csv(
  here(
    "code/analysis/12_cross-study_enrichment",
    "PEC_all_cell_types_xenium_overlap.csv"
  )
)




# Identify L6 related marker genes ----
L6_marker_df <- c("L6.CT", "L6.IT", "L6.IT.Car3", "L6b") |>
  set_names() |>
  map(
    ~ cell_type_enrich_df[
      # Rows
      which(tstats[, .x] > 0 & fdrs[, .x] < 0.05),
      # Columns
      c("ensembl", "gene")
    ] |> mutate(
      PEC_cell_type = .x
    )
  ) |>
  bind_rows()

rownames(L6_marker_df) <- NULL

# Find overlap with Xenium data ----
inner_join(L6_marker_df,
  genes_xenium,
  by = c("ensembl" = "Ensembl ID")
) |>
  filter(ensembl %in% genes_xenium$`Ensembl ID`) |>
  select(ensembl, Gene, PEC_cell_type, cell_type_updated, layer_marker, Note, everything()) |>
  select(-gene,
    "Anno_layer_marker" = "layer_marker",
    "Anno_cell_type" = "cell_type_updated",
    "Anno_cell_type_obselete" = "cell_type"
  ) |>
  # str()
  write_csv(
    here(
      "code/analysis/12_cross-study_enrichment",
      "PEC_L6_marker_genes_xenium_overlap.csv"
    )
  )
