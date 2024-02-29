# Job Script Related Notes ------------------------------------------------
# Running time: 1 hour
# Mem usage: asked for 40G, used 32GB in interactive session.



# Load Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(spatialLIBD)
  library(sessioninfo)
})


# TODO: load in spe object
raw_spe <- readRDS(
  here::here("processed-data", "rds", "01_build_spe",
             # TODO: rename
             "test_raw_spe_w_spg_N63.rds")
)

spe <- raw_spe[, raw_spe$in_tissue == 1 ]



# Visualize the distribution per sample -----------------------------------

spg_col_names <- colData(spe) |> colnames() |>
  grep("^spg_", x = _, value = TRUE)

spg_col_names

spe$spg_PBW0_5 <- as.character(spe$spg_PBW > 0.5)
spe$spg_PBW0_25 <- as.character(spe$spg_PBW > 0.25)

raw_spe$spg_PBW0_25 <- as.character(raw_spe$spg_PBW > 0.25)


saveRDS(
  data.frame(
    key = raw_spe$key,
    spg_PBW0_25 = raw_spe$spg_PBW0_25
  ),
  here("processed-data/rds/pnn/df_spot_calling.rds")
)

c("spg_NBW", "spg_PBW", "spg_CNBW" ) |> 
  lapply(
    FUN = function(.x)
      vis_grid_gene(
      spe,
      geneid = .x,
      spatial = FALSE,
      point_size = 0.8,
      pdf_file = here(
        "plots/pnn",
        paste0("test_", .x, ".pdf")
      )
    )
  )

vis_grid_clus(
  spe, clustervar  = "spg_PBW0_5",
  spatial = FALSE,
  point_size = 0.8,
  pdf_file = here(
    "plots/pnn",
    paste0("test_PBW0_5.pdf")
  )
)

vis_grid_clus(
  spe, clustervar  = "spg_PBW0_25",
  spatial = FALSE,
  point_size = 0.8,
  pdf_file = here(
    "plots/pnn",
    paste0("test_PBW25.pdf")
  )
)


print("Finish plotting spot plots")


# Plot distributions ------------------------------------------------------
## Not that helpful so far due to the excess 0s.

c("spg_NBW", "spg_PBW", "spg_CNBW" ) |> 
  lapply(
    FUN = function(.x)
      vis_grid_gene(
        spe,
        geneid = .x,
        spatial = FALSE,
        point_size = 0.8,
        pdf_file = here(
          "plots/pnn",
          paste0("test_", .x, ".pdf")
        )
      )
  )

.x <- [3]
# TODO: create the new measure that AD project used.....

ggplot() +
  geom_violin(aes(x = spe$sample_id, y= spe[[.x]])) +
  # scale_y_log10() + # O counts are removed
  labs(
    y = .x
  )







png(here("plots/pnn/pairwise_plot.png"))
seg_met_mat <- colData(spe) |> data.frame() |>
  select(c("spg_NBW", "spg_PBW", "spg_CNBW" )) |> 
  drop_na()

# Pairwise bi-variate plot ---
pairs(seg_met_mat)
dev.off()

# PCA analysis ----
pca_mdl <- prcomp(seg_met_mat, center = TRUE, scale. = TRUE)

summary(pca_mdl)

# Plotting variance explained, First two PC's good enough
plot(pca_mdl) 
# biplot(pca_mdl) # Failed because including labels


png(here("plots/pnn/test_pca_plot.png"))
ggplot(data.frame(pca_mdl$x)) +
  geom_point(aes(x = PC1, y = PC2))
dev.off()

# Thresholding Result---

## Way 1 -----------------------------------------------------------------


## Way 2 ----------------------------------------------------------------



# Session Info ------------------------------------------------------------
session_info()