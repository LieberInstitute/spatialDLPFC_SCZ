# Load library -----
suppressPackageStartupMessages({
  # library(SingleCellExperiment)
  # library(SpatialExperiment)
  library(tidyverse)
  # library(spatialLIBD)
  # library(scater)
  library(here)
})

# Load SPD data ----
## Load Spe ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)
## Load SpD data frame ----
df_spe <- readRDS(
  here(
    "processed-data", "rds", "spatial_cluster",
    "PRECAST", "test_RRECAST_label_df.rds"
  )
) |> mutate_if(is.numeric, ~ paste0("SpD_", .x))

## Merge two data frame ----

df_spe <- df_spe |>
  left_join(
    colData(spe) |> data.frame() |> select("key", "sample_id", "dx"),
    by = "key"
  )

# Calculate composition df ----

spd_var <- "PRECAST_8" # TODO: change to prefered 
df_spe$fnl_spd <- df_spe[[spd_var]]

n_spots_per_spd <- df_spe |> 
group_by(
  sample_id,
  fnl_spd
) |> 
  summarize(n_spots = n()) 


## Porportion of spots plot ----
p <- n_spots_per_spd |>
  ggplot(aes(x = sample_id, y = n_spots, fill = fnl_spd)) +
  geom_bar(position = "fill", stat = "identity") +
    guides(x = guide_axis(angle = 90)) 

p;

## Number of spots plot ----
n_spots_per_spd |>
  ggplot(aes(x = sample_id, y = n_spots, fill = fnl_spd)) +
  geom_col() +
    guides(x = guide_axis(angle = 90)) 

# proportion df -----
prop_per_spd <- n_spots_per_spd |>
  mutate(freq = n_spots / sum(n_spots))


# Session Info ----
sessioninfo::session_info()