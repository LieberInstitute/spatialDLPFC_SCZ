# Load Packages ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(tidyverse)
  library(Polychrome)
  library(sessioninfo)
})


# Load data ----
## Load spe ----
# TODO: repalce with updated version
# raw_spe <- readRDS(
#   here(
#     "processed-data/rds/01_build_spe",
#     "raw_spe_wo_SPG_N63.rds"
#   )
# )

# ## Load PRECAST df ----
# PRECAST_df <- readRDS(
#   here(
#     "processed-data/rds/spatial_cluster",
#     "PRECAST",
#     "test_clus_label_df_semi_inform_k_2-16.rds"
#   )
# )

# ## Merge PRECAST df ----
# precast_vars <- grep(
#   "^PRECAST_", colnames(PRECAST_df),
#   value = TRUE
# )
# spe <- raw_spe[, raw_spe$key %in% PRECAST_df$key]
# # raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
# col_data_df <- PRECAST_df |>
#   right_join(
#     colData(spe) |> data.frame(),
#     by = c("key"),
#     relationship = "one-to-one"
#   )
# rownames(col_data_df) <- colnames(spe)
# colData(spe) <- DataFrame(col_data_df)

spe <- readRDS(
  here(
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

precast_vars <- grep(
  "^PRECAST_", colnames(colData(spe)),
  value = TRUE
)

max_k <- 16

# Plot Clustering Results ----
precast_vars |>
  walk(
    .f = function(.x) {
      cat("Start plotting for ", .x, "\n")
      k <- unique(spe[[.x]]) |> length()

      vis_grid_clus(
        spe,
        clustervar = .x,
        sort_clust = FALSE,
        spatial = FALSE,
        pdf_file = here(
          "plots/spatial_cluster/PRECAST",
          paste0("test_semi_supervised_", .x, ".pdf")
        ),
        colors = set_names(
          Polychrome::palette36.colors(max_k)[seq.int(k)],
          unique(spe[[.x]]) |> sort()
        ),
        sample_order = unique(spe$sample_id) |> sort(),
        height = 1056,
        width = 816,
        point_size = 0.8,
        alpha = 1
      )
      cat(.x, " finished\n")
    }
  )

## Plot Legend ----
# Create a data frame with 16 categories
df <- data.frame(
  all_spd = sprintf("SpD%02d", 1:max_k)
)

# Create a vector of 16 different colors
colors <- set_names(
  Polychrome::palette36.colors(max_k),
  sprintf("SpD%02d", 1:max_k)
)

# Create the plot
tmp_p <- ggplot(df, aes(x = Category, fill = Category)) +
  geom_bar() +
  scale_fill_manual(values = colors)

library(ggpubr)
leg <- get_legend(tmp_p)

# Convert to a ggplot and print
pdf(
  here(
    "plots/spatial_cluster/PRECAST",
    paste0("test_semi_supervised_legend.pdf")
  )
)
as_ggplot(leg)
dev.off()


# Session info -----
sessioninfo::session_info()
