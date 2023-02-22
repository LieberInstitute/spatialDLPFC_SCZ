# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(cellranger) # Extract excel file
library(glue)

# Load Helper Function ----------------------------------------------------
source(
  grep(
    pattern = "/fun[^/]*\\.R",
    x = list.files(
      here("code", "raw_data_process"),
      full.names = TRUE
    ),
    value = TRUE
  )
)

# VistoSeg (split images) & Loupe --------------------------------------------------------
# Manually process
# Pre-requisite to run successfully space ranger


# Copy Raw Data File ------------------------------------------------------
# https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues/10
here("code", "raw_data_process",
     "01-create_data_folders.R") |> 
  source()

# Run Space Ranger --------------------------------------------------------
# https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues/12
here("code", "raw_data_process",
     "02-prep_for_space_ranger.R") |> 
  source()


# Nuclei Segmentation with VistoSeg ---------------------------------------
# https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues/14
# TODO: Not implemented yet


