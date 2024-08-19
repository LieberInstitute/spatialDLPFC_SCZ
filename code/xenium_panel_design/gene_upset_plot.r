library(readxl)
library(tidyverse)
library(here)
library(UpSetR)


# Load Data ----
## DEGs ----
dx_degs <- read_csv(
  here("processed-data/PB_dx_genes/test_PRECAST_07.csv")
) |>
  filter(fdr_scz <= 0.10) |>
  pull(gene)


## Layer Marker genes ----
layer_genes <- read_csv(
  here(
    "code/xenium_panel_design/",
    "layer_gene_86.csv"
  )
) |> pull(gene)


## cell_type ----
cell_type_genes <- read_xlsx(
  here(
    "code/xenium_panel_design/",
    "Xenium_SCZ_ProbeSelection_1.xlsx"
  )
) |> unlist()
cell_type_genes <- cell_type_genes[which(!is.na(cell_type_genes))]


listInput <- list(dx_degs = dx_degs, layer_genes = layer_genes, cell_type = cell_type_genes)



upset(fromList(listInput), order.by = "freq")


listInput |> unlist() |> unique() |> length()
dx_degs |> filter()
