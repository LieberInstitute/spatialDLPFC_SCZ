suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
  library(viridis)
})



# Load Dx_DEG data ---
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
)

HK_genes <- c("GUSB", "PPIA", "UBC", "SDHA", "ACTB", "TBP", "B2M", "GAPDH")

gene_df |> filter(
  gene %in% HK_genes
) |> select(
  gene, p_value_scz, fdr_scz)