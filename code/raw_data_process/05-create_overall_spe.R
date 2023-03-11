library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)


# TODO: replace this part with the preprocess_meta_file

sample_info <- data.frame(
  sample_id = list.files(
    here("processed-data", "spaceranger"),
    pattern = "^Br"
  )
) |> 
  mutate(
    sample_path = file.path(
      here::here("processed-data", "spaceranger"),
      sample_id,
      "outs"
    )
  )

stopifnot(all(file.exists(sample_info$sample_path)))

## Build SPE object
Sys.time()
spe <- spatialLIBD::read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()
