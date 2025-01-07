library(SpatialExperiment)
library(spatialLIBD)
library(here)



clustr_method <- "PRECAST"  # TODO: more options
clustr_opt <- "semi_inform_K8" # TODO: more option

spe_path <- here(
  "processed-data/rds/spatial_cluster/",
  clustr_method,
  paste0("test_PRECASTObj_", clustr_opt, ".rds")
)


spe <- readRDS(spe_path)


colnames(colData(spe))

fld_plot_spatial_cluster <- here("plots/spatial_cluster")

dir.create(
  fld_plot_spatial_cluster,
  recursive = TRUE,
  showWarnings = FALSE)



vis_grid_clus(
  spe,
  clustervar = "PRECAST_cluster",
  pdf_file = file.path(fld_plot_spatial_cluster, 
                       "PRECASTObj_semi_inform_K8.pdf")
)

