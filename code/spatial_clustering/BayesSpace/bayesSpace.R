# set up arrayjob to run k=2 to k = 15
# don't use spatial preprocess. in order to do this you have to reset metadata
# increase nrep for spatialCluster??

library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library(BayesSpace)
library("RColorBrewer")
library(ggplot2)


method_batch_correct <- "harmony"
# method_batch_correct <- "MNN"
fld_batch_correct <- here("processed-data", "rds", "spe", "batch_corrected")
spe <- readRDS(here::here(fld_batch_correct,
                          paste0("test_spe_", method_batch_correct, ".rds")))



# bayesSpace with Batch Corrected Data----

## do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
# TODO: rename the row and col to show it is concatenated.
spe$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- colData(spe)$array_col

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)
# spe_bs <- spatialPreprocess(spe, n.PCs=7,
#                             assay
#                             n.HVGs=2000, log.normalize=FALSE)


set.seed(030122)

# Make sure which 





message("Running spatialCluster()")



# TODO: adjust this
k <- 9
# k <- as.numeric(Sys.getenv("SGE_TASK_ID"))



# Pre-validation of Slices Placement --------------------------------------
fld_plot_spatialcluster <- here::here("plots", "spatial_cluster")
dir.create(file.path(fld_plot_spatialcluster, "BayesSpace"),
           showWarnings = FALSE, recursive = TRUE)



pdf(file = here::here(file.path(fld_plot_spatialcluster, "BayesSpace"),
                      "test_BayesSpace_offset_check.pdf"))
clusterPlot(spe, "sample_id", color = NA) + # make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

### BayesSpace on Batch Corrected
# TODO: make sure use the correct use.dimred for MNN
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)

spe$bayesSpace_temp <- spe$spatial.cluster
# TODO, edit the name
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name


fld_data_spatialcluster <- here("process-data", "rds", "spatial_cluster")
dir.create(
  file.path(fld_data_spatialcluster, "BayesSpace"),
  showWarnings = FALSE,
  recursive = TRUE)

saveRDS(
  colData(spe),
  file.path(
    fld_data_spatialcluster, "BayesSpace",
            paste0("test_BayesSpace_", method_batch_correct, 
                   "_", k,".rds")
            )
)


# cluster_export(
#   spe,
#   bayesSpace_name,
#   cluster_dir = here::here("processed-data", "rds", "spe", "clustering_results")
# )

# sample_ids <- unique(colData(spe)$sample_id)

# pdf(file = here::here("plots", paste0("test_vis_clus_bayesSpace_harmony_", k, ".pdf")))
# for (i in seq_along(sample_ids)) {
#   my_plot <- vis_clus(
#     spe = spe,
#     clustervar = bayesSpace_name,
#     sampleid = sample_ids[i],
#     colors = setNames(Polychrome::palette36.colors(k), seq_len(k))
#   )
#   print(my_plot)
# }
# dev.off()

### BayesSpace on Non-Batch Corrected
# spe <- spatialCluster(spe, use.dimred = "PCA", q = k, nrep = 10000)
# 
# spe$bayesSpace_temp <- spe$spatial.cluster
# bayesSpace_name <- paste0("bayesSpace_harmony_", k)
# colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name
# 
# cluster_export(
#   spe,
#   bayesSpace_name,
#   cluster_dir = here::here("processed-data", "rds", "spe", "clustering_results"),
#   overwrite = TRUE
# )
# 
# sample_ids <- unique(colData(spe)$sample_id)
# mycolors <- brewer.pal(7, "Dark2")
# 
# pdf(file = here::here("plots", paste0("test_vis_clus_bayesSpace_pca_", k, ".pdf")))
# for (i in seq_along(sample_ids)) {
#   my_plot <- vis_clus(
#     spe = spe,
#     clustervar = bayesSpace_name,
#     sampleid = sample_ids[i],
#     colors = mycolors
#   )
#   print(my_plot)
# }
# dev.off()
# 
# saveRDS(spe, here("processed-data", "rds", "spe", "test_spe_bayesSpace_harmony.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()