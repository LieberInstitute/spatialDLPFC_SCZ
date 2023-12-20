
# Load packages -----------------------------------------------------------
library(SpatialExperiment)
library(tidyverse)
library(here)


# Load raw spe object -----------------------------------------------------

spe <- readRDS(
  here("processed-data/rds/01_build_spe/",
       "raw_spe_wo_SPG_N63.rds")
)

col_df <- colData(spe) |> data.frame()


# Per-sample statistics ---------------------------------------------------
spot_df <- col_df |> 
  group_by(sample_id) |> 
  summarize(
    n_spots_in = sum(in_tissue),
    n_spots_out = sum(!in_tissue)
    ) |> 
  left_join(
    metadata(spe)$dx_df |> select(
      sample_id, dx
    ),
    by = c("sample_id" = "sample_id")
  )

# Test Statistics ---------------------------------------------------------
t.test(n_spots_in ~ dx, spot_df)
wilcox.test(n_spots_in ~ dx, spot_df)

spot_df |> 
  ggplot() +
  geom_boxplot(aes(x = dx, y = n_spots_in, color = dx))


