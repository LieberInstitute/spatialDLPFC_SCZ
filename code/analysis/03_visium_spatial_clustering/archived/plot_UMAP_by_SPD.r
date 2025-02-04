# Load library -----
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(scater)
  library(here)
})

# Load SPD data ----
## Load PRECAST SPD labels ----
### One time: save PRECAST spatial clustering res ----
# path_PRECAST_int_spe <- file.path(
#   "processed-data", "rds", "spatial_cluster",
#   "PRECAST", "test_spe_semi_inform.rds"
# )

# spe <- readRDS(path_PRECAST_int_spe)

# df_spe <- colData(spe) |>
#   data.frame() |>
#   select(key, starts_with("PRECAST"))

# saveRDS(df_spe,
#   file = here(
#     "processed-data", "rds", "spatial_cluster",
#     "PRECAST", "test_RRECAST_label_df.rds"
#   )
# )

df_spe <- readRDS(
  here(
    "processed-data", "rds", "spatial_cluster",
    "PRECAST", "test_RRECAST_label_df.rds"
  )
) |> mutate_if(is.numeric, ~ paste0("SpD_", .x))




## Load UMAP data ----
spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "test_raw_spe_UMAP_N63.rds"
  )
)

## Merge SPD labels with UMAP spe ----
df_merg <- colData(spe) |>
  data.frame() |>
  left_join(
    y = df_spe,
    by = "key"
  )

rownames(df_merg) <- df_merg$key

colData(spe) <- df_merg |> DataFrame()

# Plot UMAP with SpD ----

var_spd <- "PRECAST_7"
.dimred <- "UMAP_2000"
fld_path <- "plots/spatial_cluster/UMAP"

stopifnot(.dimred %in% reducedDimNames(spe))
stopifnot(file.exists(fld_path))

## Combined Plot ----
png(here(
  fld_path,
  paste0(.dimred, "_", var_spd, ".png")
))
plotReducedDim(spe,
  dimred = .dimred,
  point_size = 0.1,
  colour_by = var_spd
) +
  theme(legend.position = "none") #+
# scale_colour_manual(values = colors)
dev.off()

# Sample Specific UMAP
# for (i in 1:3){
for (i in seq.int(unique(spe[[var_spd]]))) {
  .smp <- unique(spe[[var_spd]])[i]
  spe$sub_sample <- FALSE
  spe$sub_sample[spe[[var_spd]] == .smp] <- TRUE

  # order(spe$sub_sample)

  png(
    here(
      fld_path,
      paste0(.dimred, "_", var_spd, "_", .smp, ".png")
    )
  )
  print(
    plotReducedDim(
      spe[, order(spe$sub_sample)],
      dimred = .dimred,
      point_size = 0.1,
      colour_by = "sub_sample"
    ) +
      theme(legend.position = "none") +
      scale_colour_manual(
        values = c("FALSE" = "lightgrey", "TRUE" = "red")
      ) +
      labs(title = .smp)
  )
  dev.off()
}
