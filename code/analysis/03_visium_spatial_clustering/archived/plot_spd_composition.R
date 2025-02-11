# Load library -----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(tidyverse)
  # library(spatialLIBD)
  # library(scater)
  library(here)
  library(sessioninfo)
})

# Load SPD data ----
## Load Spe ----
spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# ## Load SpD data frame ----
# df_spe <- readRDS(
#   here(
#     "processed-data", "rds", "spatial_cluster",
#     "PRECAST", "test_RRECAST_label_df.rds"
#   )
# ) |> mutate_if(is.numeric, ~ paste0("SpD_", .x))

# ## Merge two data frame ----

df_spe <- spe |>
  colData() |>
  data.frame()

# Calculate composition df ----

spd_var <- "PRECAST_07" # TODO: change to prefered
df_spe$fnl_spd <- df_spe[[spd_var]]

n_spots_per_spd <- df_spe |>
  group_by(
    sample_id,
    fnl_spd
  ) |>
  summarize(n_spots = n())


n_spots_per_spd <- n_spots_per_spd |> left_join(
  metadata(spe)$dx_df |> select(sample_id, dx, subject),
  by = "sample_id"
)

## Porportion of spots plot ----
pdf(
  here("plots", "spatial_cluster",
  "spd_comp_PRECAST_07.pdf"
  )
)
p <- n_spots_per_spd |>
  ggplot(aes(x = subject, y = n_spots, fill = fnl_spd)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~dx, scales = "free_x") +
  scale_fill_manual(
    name = "Spatial Domain",
    values = set_names(
          Polychrome::palette36.colors(16)[seq.int(7)],
          unique(n_spots_per_spd$fnl_spd) |> sort()
        )
  ) +
  labs(
    y = "Prop. of spots"
  ) +
  guides(
    x = guide_axis(angle = 90)
  ) +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 20)
  )
print(p)
dev.off()

## Number of spots plot ----
# n_spots_per_spd |>
#   ggplot(aes(x = sample_id, y = n_spots, fill = fnl_spd)) +
#   geom_col() +
#   guides(x = guide_axis(angle = 90))



# Session Info ----
sessioninfo::session_info()
