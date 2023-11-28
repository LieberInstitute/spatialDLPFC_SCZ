# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(gtsummary)


path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

spe <- readRDS(
  path_spe_after_spot_qc
)

col_df <- colData(spe) |> data.frame()


# Tables ------------------------------------------------------------------
sample_mean_df <- col_df |> 
  select(sample_id, 
         sum_umi, sum_gene, expr_chrM, 
         expr_chrM_ratio) |> 
  group_by(sample_id) |> 
  summarize(
    # Mean stat
    across(sum_umi:expr_chrM_ratio, mean),
    # Number of spots
    n = n()
  ) |> 
  ungroup()

demo_df <- metadata(spe)$dx_df

sample_mean_df <- sample_mean_df |> 
  left_join(
    y = demo_df |> select(sample_name, dx),
    by = c("sample_id" = "sample_name")
  )

# TODO: why there are NaN value
stopifnot("Missing stat exists, please check" =
            sum(!complete.cases(sample_mean_df)) == 0)
# Investigation code
# sample_mean_df |> filter(!complete.cases(sample_mean_df))


# * Per-sample ------------------------------------------------------------
# TODO: Needs formating
col_df |> 
  select(sample_id, 
         sum_umi, sum_gene, expr_chrM, 
         expr_chrM_ratio) |> 
  tbl_summary(
    by = sample_id
  )

# * Per-dx ----------------------------------------------------------------
sample_mean_df |> 
  select(-sample_id) |> 
  tbl_summary(
    by = dx
  ) |> 
  add_p()


# Plot --------------------------------------------------------------------
ggplot(col_df|> 
         left_join(
           y = demo_df |> select(sample_name, dx),
           by = c("sample_id" = "sample_name")
         )) + 
  geom_violin(aes(x = sample_id, y = sum_umi)) +
  facet_grid(. ~ dx,  scales = "free_x") +
  scale_y_log10() +
  scale_x_discrete(guide = guide_axis(angle = 90))

# TODO:
# 1) save plot
# 2) add some sort of p-values if necessary.




