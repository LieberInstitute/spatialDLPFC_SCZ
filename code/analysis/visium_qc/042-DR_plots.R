# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(scater)
library(scran)
# library(pryr)                 # Check spe size
# library(spatialLIBD)
library(tidyverse)


# TODO: edit
path_spe_after_spot_qc <- here::here("processed-data", "rds", "spe", 
                                     "spe_after_spot_qc.rds")


spe <- readRDS(
  path_spe_after_spot_qc
)

# Create logcounts
spe <- logNormCounts(spe)

spe$slide_num <- str_remove(spe$sample_id, "_[A-D][1]")

# For older version of metadata(spe)
# spe$dx <- metadata(spe)$dx[match(spe$sample_id, metadata(spe)$sample_name)]
spe$dx <- metadata(spe)$dx_df$dx[match(spe$sample_id, metadata(spe)$sample_name)]


# PCA ---------------------------------------------------------------------
spe <- runPCA(spe) # TODO: edit this
plotPCA(spe, color_by = "sample_id")
plotPCA(spe, color_by = "slide_num")
plotReducedDim(spe, dimred="PCA", ncomponents=4,
               color_by = "slide_num")

percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")


# Default Walktrap method
# n=50
walktrap_50 <- clusterCells(spe, use.dimred="PCA")
table(walktrap_50)
spe$walktrap_PCA_50 <- walktrap_50

colLabels(spe) <- nn.clusters

plotPCA(spe, color_by = "walktrap_PCA_50")

table(spe$label, spe$sample_id)

# PCA n=10
spe <- runPCA(spe, ncomponents = 10, name = "PCA_10")
walktrap_10 <- clusterCells(spe, use.dimred="PCA_10")
spe$walktrap_PCA_10 <- walktrap_10
plotReducedDim(spe, dimred="PCA_10", ncomponents=2, color_by = "walktrap_PCA_10")


library(spatialLIBD)
# Bug when plotting V12F14-057_A1
vis_grid_clus(spe,
              clustervar = "walktrap_PCA_50",
              pdf_file = here("plots/01_build_spe/test_walk_trap_50.pdf")
)

vis_grid_clus(
  spe[, spe$sample_id =="V12F14-057_A1"],
  clustervar = "walktrap_PCA_50",
  return_plots = TRUE
  # pdf_file = here("plots/01_build_spe/test_walk_trap_50.pdf")
)

vis_grid_clus(spe,
              clustervar = "walktrap_PCA_10",
              pdf_file = here("plots/01_build_spe/test_walk_trap_10.pdf")
)



# plotReducedDim(spe, dimred="PCA", ncomponents=4,
#                color_by = "label")
# 
# plotReducedDim(spe, dimred="PCA", ncomponents=4,
#                color_by = "slide_num")


# TSNE --------------------------------------------------------------------
set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=20)





# UMAP --------------------------------------------------------------------
set.seed(2)
spe <- runUMAP(spe, dimred = "PCA") # TODO: edit this
# Same as below
# # set.seed(1)
# spe <- runUMAP(spe, dimred = "PCA", n_dimred = 10, name = "UMAP_10PCs_auto")
set.seed(1)
spe <- runUMAP(spe, dimred = "PCA_10", n_dimred = 10, name = "UMAP_10PCs")# TODO: edit this
spe <- runUMAP(spe, dimred = "PCA_10", n_dimred = 5, name = "UMAP_5PCs")# TODO: edit this

# ggpubr::ggarrange(
#   plotReducedDim(spe, dimred="UMAP_10PCs_auto", color_by = "sample_id"),
#   plotReducedDim(spe, dimred="UMAP_10PCs", color_by = "sample_id")
# )



ggpubr::ggarrange(
  plotUMAP(spe, color_by = "sample_id"), 
  plotReducedDim(spe, dimred="UMAP_10PCs", color_by = "sample_id"),
  common.legend = TRUE,
  legend = "bottom"
)

ggpubr::ggarrange(
  plotUMAP(spe, color_by = "sample_id"), 
  plotUMAP(spe, color_by = "walktrap_PCA_50"),
  # common.legend = TRUE,
  legend = "bottom"
)

plotReducedDim(spe, dimred="UMAP_5PCs", color_by = "sample_id")


plotUMAP(spe[, spe$sample_id==c("V12F14-057_A1", "V12F14-057_C1")], color_by = "sample_id")
plotUMAP(spe, color_by = "slide_num")
plotUMAP(spe, color_by = "dx")
