# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)


path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

spe <- readRDS(
  path_spe_after_spot_qc
)


# Range of SPG channels ---------------------------------------------------

col_df <- colData(spe) |> data.frame()

col_df |> group_by(sample_id) |> 
  summarize(
    across(starts_with("spg_N"),
           list (#"min" = min,  # Seems all mins are 0 
                 "max" = max))
  )

# sample_id     spg_NAF_max spg_NClaudin5_max spg_NDAPI_max spg_NNeuN_max spg_NWFA_max
# <glue>              <int>             <int>         <int>         <int>        <int>
# 1 V12F14-053_A1          47                18            20            13           41
# 2 V12F14-053_B1          71                25            18            18           19
# 3 V12F14-053_C1          69                20            27            16           37
# 4 V12F14-053_D1          50                29            20            77           40
# 5 V12F14-057_A1         253                32            18            24           12
# 6 V12F14-057_B1          56                24            24            30           44
# 7 V12F14-057_C1          88                33            26           112           34
# 8 V12F14-057_D1          71                21            17            33           64



#  Per-sample Plots -----------------------------------------------

names(col_df) |> 
  grep("^spg_N", x = _,
       value = TRUE) |> 
  walk(
    ~ vis_grid_gene(
      spe = spe,
      geneid = .x,
      # geneid = "spg_NDAPI",
      pdf = here::here("plots", "02_visium_qc", "SPG",
                       paste0("spot_plot_", 
                              str_remove(.x, "spg_N"), ".pdf")),
      spatial=FALSE,
      assayname = "counts"
    )
  )
  

