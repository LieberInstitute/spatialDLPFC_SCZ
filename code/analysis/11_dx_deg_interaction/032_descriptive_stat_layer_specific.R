# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load data ----
## Load Layer-ajusted DEGs ----
spd_deg_df <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_restricted_degs_all_spds.csv"
  )
)


# Descriptive Stats ----


## Layer-restricted DEGs ----
spd_deg_df |>
  filter(P.Value < 0.05) |>
  summarize(n = n())

spd_deg_df |>
  filter(P.Value < 0.05) |>
  summarize(n = n_distinct(gene_id))


spd_deg_df |>
  filter(adj.P.Val < 0.10) |>
  summarize(n = n())

spd_deg_df |>
  filter(adj.P.Val < 0.10) |>
  summarize(n = n_distinct(gene_id))

# NOTE: also see the file `code/analysis/11_dx_deg_interaction/031_descriptive_stat_layer_restricted_genes.R`


## Layer-specific DEGs ----
### Total number ----
spd_deg_df |>
  filter(layer_specific) |>
  summarize(n = n())

### Total number by up-/down-reg ----
spd_deg_df |>
  filter(layer_specific) |>
  mutate(logFC_sign = if_else(logFC > 0, "up", "down")) |>
  group_by(logFC_sign) |>
  summarize(n = n())


### PRECAST_spd ----
spd_deg_df |>
  filter(layer_specific) |>
  group_by(PRECAST_spd) |>
  summarize(n = n())

### Call out indiviual genes in layer-restricted DEGs ----

spd_deg_df |>
  filter(layer_specific) |>
  filter(gene %in% c("NRG4"))


spd_deg_df |>
  filter(layer_specific) |>
  filter(gene %in% c(
    "KCNE4", "PDLIM5", "GNG12", "EEF2K",
    "NR4A1", "RPL26", "SMAD6", "TBX3",
    "LILRB1", "HAVCR2", "CCL2", "HLA-DMB", "HLA-DQB1", "CD44"
  ))
