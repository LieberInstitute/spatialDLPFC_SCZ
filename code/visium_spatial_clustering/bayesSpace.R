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

spe <- readRDS(here::here("processed-data", "rds",
                          "spe", "01_build_spe",
                          "test_spe_filtered_final.Rdata"))

# TODO: adjust this
k <- 9
# k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(030122)

# TODO: change path
pdf(file = here::here("plots",  "test_BayesSpace_offset_check.pdf"))
clusterPlot(spe, "sample_id", color = NA) + # make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

### BayesSpace on Batch Corrected
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rds", "spe", "clustering_results")
)

sample_ids <- unique(colData(spe)$sample_id)

pdf(file = here::here("plots", paste0("test_vis_clus_bayesSpace_harmony_", k, ".pdf")))
for (i in seq_along(sample_ids)) {
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors = setNames(Polychrome::palette36.colors(k), seq_len(k))
  )
  print(my_plot)
}
dev.off()

### BayesSpace on Non-Batch Corrected
spe <- spatialCluster(spe, use.dimred = "PCA", q = k, nrep = 10000)

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rds", "spe", "clustering_results"),
  overwrite = TRUE
)

sample_ids <- unique(colData(spe)$sample_id)
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots", paste0("test_vis_clus_bayesSpace_pca_", k, ".pdf")))
for (i in seq_along(sample_ids)) {
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors = mycolors
  )
  print(my_plot)
}
dev.off()

saveRDS(spe, here("processed-data", "rds", "spe", "test_spe_bayesSpace_harmony.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()