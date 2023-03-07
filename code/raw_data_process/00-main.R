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

# Space Ranger QC ---------------------------------------------------------
# https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues/16

# TODO:  Add here


# Nuclei Segmentation with VistoSeg ---------------------------------------
# https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues/14
here("code", "raw_data_process",
     "04-vistoseg.R") |> 
  source()



