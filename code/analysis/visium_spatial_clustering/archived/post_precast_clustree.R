# Load Libray -------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  # library(spatialLIBD)
  library(tidyverse)
  # library(escheR)
  # library(ggpubr)
  library(mclust)
  library(clustree)
  library(sessioninfo)
})


# Path --------------------------------------------------------------------
fld_data_spatialcluster <- here(
  "processed-data",
  "rds", "spatial_cluster")

path_PRECAST_int_spe <- file.path(
  fld_data_spatialcluster, "PRECAST",
  paste0("test_spe_semi_inform",".rds")
)

# Load data ---------------------------------------------------------------
spe <- readRDS(
  path_PRECAST_int_spe
)

spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )]

pdf(here("plots/spatial_cluster/assessment/test_clustree.pdf"))
# Overall -----------------------------------------------------------------
(
  clustree(colData(spe) |> data.frame(), prefix = "PRECAST_") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +
  theme(legend.position = "bottom") +
  labs(title = "All Sample")
) |> print()


# Per Dx Group ----
for(.dx in unique(spe$dx)){
  spe_sub <- spe[, spe$dx == .dx]
  (
    clustree(colData(spe_sub) |> data.frame(), prefix = "PRECAST_") +
      guides(edge_colour = FALSE, edge_alpha = FALSE) +
      theme(legend.position = "bottom") +
      labs(title = paste0("spe$dx = ",.dx))
  ) |> print()
}




# Sample Specific ---------------------------------------------------------

for(.smp_name in unique(spe$sample_id)){
  spe_sub <- spe[, spe$sample_id == .smp_name]
  (
    clustree(colData(spe_sub) |> data.frame(), prefix = "PRECAST_") +
      guides(edge_colour = FALSE, edge_alpha = FALSE) +
      theme(legend.position = "bottom") +
      labs(title = .smp_name)
  ) |> print()
}



dev.off()



# SessionInfo -------------------------------------------------------------
sessioninfo::session_info()

