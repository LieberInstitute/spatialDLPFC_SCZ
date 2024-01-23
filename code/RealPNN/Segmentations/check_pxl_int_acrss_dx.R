library(tidyverse)
library(gtsummary)
library(here)

df <- read_csv(here("processed-data/VistoSeg/captureAreas/image_metrics.csv"))


# TODO: read in sample

dx_df <- metadata(spe)$dx_df

full_df <- df |> left_join(
  y = dx_df,
  by = c("Sample" = "sample_id")
)

full_df |> 
  tbl_summary(
    include = c(Mean, Median, Mode, Kurtosis, Width ),
    by = dx) |> 
  add_n()  |>  # add column with total number of non-missing observations
  add_p() # test for a difference between groups
