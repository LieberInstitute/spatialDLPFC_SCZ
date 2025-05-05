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


# Remove first two slides
sub_df <- full_df |> 
  filter(str_starts(Sample, "V12F14-053", negate = TRUE)) |> 
  filter(str_starts(Sample, "V12F14-057", negate = TRUE))
sub_df |>  
  tbl_summary(
    include = c(Mean, Median, Kurtosis, Width ),
    by = dx) |> 
  add_n()  |>  # add column with total number of non-missing observations
  add_p() # test for a difference between groups

ggplot(sub_df) +
  geom_jitter(aes(x = dx, y = Width))

