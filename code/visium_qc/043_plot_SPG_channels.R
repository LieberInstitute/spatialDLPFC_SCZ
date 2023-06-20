# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)


path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

# In-tissue spots only
spe <- readRDS(
  path_spe_after_spot_qc
)

# TODO: remove the log normalization, should be done in prev step
# library(scran)
# spe <- logNormCounts(spe)

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
  grep("^spg_[P|N]", x = _,
       value = TRUE) |> 
  walk(
    ~ vis_grid_gene(
      spe = spe,
      geneid = .x,
      # geneid = "spg_NDAPI",
      pdf = here::here("plots", "02_visium_qc", "SPG",
                       paste0("spot_plot_", 
                              str_remove(.x, "spg_"), ".pdf")),
      spatial=FALSE,
      assayname = "counts"
    )
  )


# Plot for pvalb to confirm WFA seg ---------------------------------------


vis_grid_gene(
  spe = spe,
  geneid = "ENSG00000100362",
  # geneid = "spg_NDAPI",
  pdf = here::here("plots", "02_visium_qc", "SPG",
                   paste0("spot_plot_", 
                          "PVALB_counts", ".pdf")),
  spatial=FALSE,
  assayname = "counts"
)

# Plot for PCP4 to confirm orientation seg ---------------------------------------


vis_grid_gene(
  spe = spe,
  geneid = "ENSG00000183036",
  # geneid = "spg_NDAPI",
  pdf = here::here("plots", "02_visium_qc", "SPG",
                   paste0("spot_plot_", 
                          "PCP4_counts", ".pdf")),
  spatial=FALSE,
  assayname = "counts"
)

vis_grid_gene(
  spe = spe,
  geneid = "ENSG00000183036",
  # geneid = "spg_NDAPI",
  pdf = here::here("plots", "02_visium_qc", "SPG",
                   paste0("spot_plot_", 
                          "PCP4_logcounts", ".pdf")),
  spatial=FALSE,
  assayname = "logcounts"
)


# WFA - PVAL Quantative Investigation -------------------------------------




stopifnot(length(counts(spe)['ENSG00000100362', ]) == nrow(col_df))
col_df$PVALB_count <- counts(spe)['ENSG00000100362', ]
col_df$PVALB_logcount <- logcounts(spe)['ENSG00000100362', ]


# Sparsity
col_df |> group_by(sample_id) |> 
  summarize(
    n = n(),
    NWFA_pos = sum(spg_NWFA > 0),
    PWFA_pos = sum(spg_PWFA > 0),
    # NOTE: NWFA and PWFA is the same
    PVALB_pos = sum(PVALB_count > 0)
  ) |> 
  ungroup() |> 
  mutate(
    # Proportion of spots contains PNN (WFA) segmentation
    WFA_spars_prop = NWFA_pos/n,
    PVALB_spars_prop = PVALB_pos/n
  )

# NOTE: The sparsity should not be directly comparable (column-wise comparison)
#       due to morphology composition variation.

# sample_id         n NWFA_pos PWFA_pos PVALB_pos WFA_spars_prop PVALB_spars_prop
# <glue>        <int>    <int>    <int>     <int>          <dbl>            <dbl>
# 1 V12F14-053_A1  4933     1567     1567       930         0.318            0.189 
# 2 V12F14-053_B1  4238      996      996      1110         0.235            0.262 
# 3 V12F14-053_C1  3612     1361     1361       940         0.377            0.260 
# 4 V12F14-053_D1  3189      587      587       789         0.184            0.247 
# 5 V12F14-057_A1  3099       69       69       139         0.0223           0.0449
# 6 V12F14-057_B1  3815     1006     1006       634         0.264            0.166 
# 7 V12F14-057_C1  3288     1079     1079       499         0.328            0.152 
# 8 V12F14-057_D1  4069     1589     1589       915         0.391            0.225

# NOTE: the comparison between WFA_spars_prop and PVALB_spars_prop provides
#       indicator which sample's segmentation could be bad
# Seems like WFA_spars_prop should be larger than PVALB_spars prop.

# Tissue morphology is related to the overall segmentation performance





# Box Plot range of WFA and PWFA among spots
ggplot(col_df |> filter(spg_NWFA > 0)) +
  geom_boxplot(aes(x = sample_id, y = spg_NWFA))

# TODO: possible QC with PWFA
ggplot(col_df |> filter(spg_PWFA > 0)) +
  geom_boxplot(aes(x = sample_id, y = spg_PWFA))
  
ggplot(col_df |> filter(PVALB_count > 0)) +
  geom_boxplot(aes(x = sample_id, y = PVALB_count))

ggplot(col_df |> filter(PVALB_logcount > 0)) +
  geom_boxplot(aes(x = sample_id, y = PVALB_logcount))



# Plot range of PVAL

  


# V12F14−057_A1 -----------------------------------------------------------
channel_plot_list <- names(col_df) |> 
  grep("^spg_N", x = _,
       value = TRUE) |> 
  map(
    ~ vis_gene(
      spe = spe[, spe$sample_id=="V12F14-057_A1"],
      geneid = .x,
      # geneid = "spg_NDAPI",
      # pdf = here::here("plots", "02_visium_qc", "SPG",
      #                  paste0("spot_plot_", 
      #                         str_remove(.x, "spg_N"), ".pdf")),
      spatial=FALSE,
      assayname = "counts",
      return_plots = TRUE
    )
  )

cowplot::plot_grid(plotlist = channel_plot_list, ncol = 3) |> 
  ggsave(here::here("plots", "02_visium_qc", "SPG",
                    paste0("spot_plot_V12F14−057_A1_all_channels.pdf")), 
           plot = _
         )
