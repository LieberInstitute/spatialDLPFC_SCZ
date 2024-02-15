library(tidyverse)
library(here)


fld_samples <- list.files(
  path = here("processed-data/spaceranger"),
  pattern = "^V.*"
  )


# fld_samples

PNN_seg_path <- file.path(
  here("processed-data/spaceranger"),
  fld_samples,
  "outs/spatial/tissue_spot_counts.csv"
) |> 
  purrr::set_names(fld_samples)

PNN_seg_path |> file.exists() |> sum()

PNN_seg_path <- PNN_seg_path[PNN_seg_path |> file.exists()]

head(PNN_seg_path)
