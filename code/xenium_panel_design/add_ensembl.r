# Load Packages ----
library(tidyverse)
library(readxl)
library(here)
library(sessioninfo)


# Load Probe file ----
prob_df <- read_xlsx(
  here("code/xenium_panel_design/Xenium_SCZ_ProbeSelection5_SHK_v1.xlsx"),
  sheet = "SHK_Final_300_v1",
  col_names = FALSE
) |> select(-`...6`)
## Remove useless columns ----
colnames(prob_df) <- c("gene", "cell_type", "deg", "layer", "other")


# Read gene info ----
raw_spe <- readRDS(
  here::here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)

rowData(raw_spe) |> colnames()

prob_df$ensembl <- rowData(raw_spe)$gene_id[
  match(prob_df$gene, rowData(raw_spe)$gene_name)
]

stopifnot(sum(is.na(prob_df$ensembl))==0)

write_csv(
  prob_df,
  here("code/xenium_panel_design/Xenium_SCZ_ProbeSelection5_SHK_v1_w_ensembl.csv"),
  na = ""
)
