# Load packages----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
})

# Load Data ----
spe <- readRDS(
  here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_w_spg_N63.rds"
  )
)

# Call SPG spots ----
spe$PNN <- ifelse(spe$spg_PWFA > 0.05, "PNN+", "PNN-")
spe$DAPI <- ifelse(spe$spg_PDAPI > 0.05 & spe$spg_PDAPI < 0.5, "DAPI+", "DAPI-")
spe$NeuN <- ifelse(spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3, "NeuN+", "NeuN-")
spe$Claudin <- ifelse(
  spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,
  "Claudin+", "Claudin-"
)

# Load SpD information ----
finalized_spd <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

col_data_df <- colData(spe) |>
  data.frame() |>
  left_join(
    finalized_spd,
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)

# Make bar plots ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))


## Proportion of the PNN+ spots per layer ----
col_dat <- colData(spe) |> data.frame()

col_dat |>
  group_by(sample_id, PRECAST_07) |>
  summarize(
    prop = sum(pnn_pos) / n()
  ) |>
  left_join(
    metadata(spe)$dx_df |> select(sample_id, dx),
    by = "sample_id"
  ) |>
  ggplot() +
  # geom_point(aes(x = PRECAST_07, y = concord, color = dx, group = dx)) +
  geom_boxplot(aes(x = PRECAST_07, y = prop, color = dx)) +
  scale_x_discrete(
    limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
    labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
  ) +
  theme_minimal()

## Proportion of the PNN+ spots ----

col_dat |>
  group_by(sample_id) |>
  summarize(
    pnn_prop = sum(pnn_pos) / n()
  ) |>
  left_join(
    metadata(spe)$dx_df |> select(sample_id, dx),
    by = "sample_id"
  ) |>
  ggplot() +
  geom_boxplot(
    aes(x = dx, y = pnn_prop, color = dx)
  ) +
  geom_jitter(
    aes(x = dx, y = pnn_prop),
    size = 0.6,
    alpha = 0.4
  ) +
  theme_minimal()
