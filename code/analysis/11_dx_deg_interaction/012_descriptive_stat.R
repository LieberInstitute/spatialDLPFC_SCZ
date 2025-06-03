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
