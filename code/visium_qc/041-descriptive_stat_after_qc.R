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

col_df <- colData(spe) |> data.frame()
library(gtsummary)

col_df |> 
  select(sample_id, 
         sum_umi, sum_gene, expr_chrM, 
         expr_chrM_ratio) |> 
  tbl_summary(
    by = sample_id
  )

