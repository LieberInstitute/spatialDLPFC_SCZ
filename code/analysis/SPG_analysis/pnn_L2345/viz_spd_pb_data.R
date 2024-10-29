# Load packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
  library(SingleCellExperiment)
})

# Load Spatial domain PB Data ----
sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "PB_dx_spg",
    "test_SPD_pseudo_pnn_pos.rds"
  )
)


col_df <- colData(sce_pseudo) |> data.frame()

# Show which domains are included ----

col_df |>
  ggplot() +
  geom_point(aes(x = PRECAST_07, y = sample_id, size = ncells)) +
  theme_minimal() +
  scale_x_discrete(
    limits = c("spd07", "spd06", "spd02", "spd05", "spd03", "spd01", "spd04")
  )


# Create SpD only data ------

## L2/3 (spd06), L3/4 (spd02), and L5 (spd05) only ----
for (.spd in c("spd06", "spd02", "spd05")) {
  sce_pseudo[, sce_pseudo$registration_variable == .spd] |>
    saveRDS(here(
      "processed-data", "rds", "PB_dx_spg",
      paste0("test_donor_pseudo_pnn_pos_", .spd, ".rds")
    ))
}



##  L2/3 (spd06), L3/4 (spd02), and L5 (spd05) together ---
sce_pseudo[, sce_pseudo$registration_variable %in% c("spd06", "spd02", "spd05")] |>
  saveRDS(here(
    "processed-data", "rds", "PB_dx_spg",
    "test_SPD_pseudo_pnn_pos_L2345.rds"
  ))





# Session info ----
session_info()
