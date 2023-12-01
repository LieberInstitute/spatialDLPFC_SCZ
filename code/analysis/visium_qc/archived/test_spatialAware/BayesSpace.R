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
                          "spe_raw.rds"))

# TODO: adjust this
k <- 2
# k <- as.numeric(Sys.getenv("SGE_TASK_ID"))



sample_id <- unique(spe$sample_id)[2]

sub_spe <- spe[, spe$sample_id == sample_id]

# Create Sample specific PCA

sub_spe$row <- sub_spe$array_row
sub_spe$col <- sub_spe$array_col

# Run BayesSpace
set.seed(030122)
obj_BS <- spatialPreprocess(sub_spe,
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE)

obj_BS <- spatialCluster(obj_BS, # TODO: change this
                        q = k, nrep = 10000)

# Convert back to running SPE.

sub_spe$BayesSpace_K2 <- obj_BS$spatial.cluster


# Visualize Result

my_plot <- vis_clus(
  spe = sub_spe,
  clustervar = "BayesSpace_K2",
  sampleid = sample_id,
  colors = setNames(Polychrome::palette36.colors(k), seq_len(k))
)


