library(readxl)
library(here)
library(UpSetR)


# Load Data ----
## DEGs ----
dx_degs <- read_csv(
  here("processed-data/PB_dx_genes/test_PRECAST_07.csv")
) |> filter(fdr_scz <= 0.10) |> pull(ensembl)


## Layer Marker genes ----
gene_names


## Xenium Brain Panel ----
panel_genes <- read_xlsx(here("code/xenium_panel_design/Xenium Human Brain Panel Gene List.xlsx")) |> pull(`Ensembl ID`)

listInput <- list(dx_degs = dx_degs, layer_genes = gene_names, panel = panel_genes)



movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
    header = T, sep = ";")



upset(fromList(listInput), order.by = "freq")


intersect(dx_degs, gene_names)
dx_degs |> filter()