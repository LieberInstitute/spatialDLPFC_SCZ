# Load packages ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(tidyverse)
  library(ggbeeswarm)
  # library(escheR)
  library(sessioninfo)
})

# Load Data ----
## Load SPG spe object ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

## Call SPG + spots ----
# Call SPG spots ----
spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
# NOTE: neuropil spot are spots doesn't have DAPI staining
spe$neuropil_pos <- ifelse(
  spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5,
  FALSE, TRUE
)
spe$neun_pos <- ifelse(
  spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,
  TRUE, FALSE
)
spe$vasc_pos <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  TRUE, FALSE
)


# EDA ----
# # of spots that are both neun+ and vacs+
table(spe$neun_pos, spe$vasc_pos)
#  3478/279806
# [1] 0.01243004

# Is this number vary by donor?
tmp_df <- colData(spe) |>
  data.frame() |>
  group_by(sample_id) |>
  summarize(
    n_both = sum(neun_pos & vasc_pos),
    prop_both = n_both / n(),
    dx = unique(dx)
  ) |>
  arrange(desc(prop_both));

# Match with the previous number
tmp_df$n_both |> sum()

head(tmp_df, n=20)

tmp_df |> 
ggplot(aes(x = dx, y = prop_both, fill = dx)) +
# geom_boxplot() +
# geom_violin() +
# geom_jitter(aes(color = dx), alpha = 0.7)+
geom_quasirandom(aes(color = dx), alpha = 0.7)

# NOTE: for this plot, it looks like the prop for both are uniform for ntc and not that uniform in scz. This observation is more obvioud in a violin plot

tmp_df |> 
ggplot(aes(x = dx, y = prop_both, fill = dx)) +
# geom_boxplot() +
geom_violin()+
geom_jitter(alpha = 0.7)

# Session Info ----
sessioninfo::session_info()
