# Load library ----
suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(tidyver)
})

# Load Spd spatial experiment object ----
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)


spe <- readRDS(
  path_PRECAST_int_spe
)

# Load spot calling ----
df_PNN <- readRDS(
    here("processed-data/rds/pnn/df_spot_calling.rds")
)

spe$PNN_spot <- df_PNN$spg_PBW0_25[match(spe$key, df_PNN$key)] |>
                    factor()

# (Only test) subsetting sample that spot calling exists ----

# colnames(colData(spe))






# Heatmap ---- 

sum_df <- data.frame(
  sample_id = spe$sample_id,
  spd = paste0("spd_", spe$PRECAST_9),
  PNN_spot = spe$PNN_spot
) |> group_by(
  sample_id, spd
) |>
summarize(
  n = n(),
  prop_true = sum(PNN_spot == "TRUE") / n()
)


# TODO: maybe change this to 
png(here("plots/pnn/prop_PNN_spd_sample.png"))
print(ggplot(sum_df, aes(spd, sample_id)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = prop_true))
)
dev.off()


ggplot(sum_df) +
  geom_point(aes(x = .data$spd, y = .data$sample_id, 
        size = .data$prop_true, col = .data$n))
# Session Info ----
sessioninfo::session_info()