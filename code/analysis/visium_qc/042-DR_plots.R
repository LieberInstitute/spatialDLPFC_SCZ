
# JHPCE Interactive -------------------------------------------------------------------
# srun --pty --x11 --mem=80G bash


# Load Packages -----------------------------------------------------------Advanced SLURM/Cluster topics
library(here)
library(SpatialExperiment)
library(scater)
library(scran) 
# library(pryr)                 # Check spe size
# library(spatialLIBD)
library(tidyverse)

# Load QC-ed Data ---------------------------------------------------------
path_spe_after_spot_qc <- here::here(
  "processed-data/rds/",
  #TODO: replace this path
  "test_spe_after_spot_qc_63.rds"
)

spe <- readRDS(
  path_spe_after_spot_qc
)

# Create logcounts
if(! "logcounts" %in% assayNames(spe))
  spe <- logNormCounts(spe)


# All Genes -------------------------------------------------------
## PCA ---------------------------------------------------------------------
set.seed(1)
spe_all <- runPCA(spe) # TODO: edit this
# spe_all <- readRDS("~/spe_pca.rds")
# Size is too big. 
# TODO: go to cluster and export only PCA assay.



# plotPCA(spe_all, color_by = "sample_id") |>
# ggsave("~/spe_all_pca_sample.pdf")
# plotPCA(spe_all, color_by = "slide_num")
plotReducedDim(spe_all, dimred="PCA", ncomponents=2,
               color_by = "sample_id")

percent.var <- attr(reducedDim(spe_all, "PCA"), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")


### (Deprecated) Clustering PCA ------------------------------
# # n=50
# walktrap_50 <- clusterCells(spe, use.dimred="PCA")
# table(walktrap_50)
# spe$walktrap_PCA_50 <- walktrap_50
# colLabels(spe) <- nn.clusters
# plotPCA(spe, color_by = "walktrap_PCA_50")
# table(spe$label, spe$sample_id)
# 
# # PCA n=10
# spe <- runPCA(spe, ncomponents = 10, name = "PCA_10")
# walktrap_10 <- clusterCells(spe, use.dimred="PCA_10")
# spe$walktrap_PCA_10 <- walktrap_10
# plotReducedDim(spe, dimred="PCA_10", ncomponents=2, color_by = "walktrap_PCA_10")
# library(spatialLIBD)
# # Bug when plotting V12F14-057_A1
# vis_grid_clus(spe,
#               clustervar = "walktrap_PCA_50",
#               pdf_file = here("plots/01_build_spe/test_walk_trap_50.pdf")
# )
# 
# vis_grid_clus(
#   spe[, spe$sample_id =="V12F14-057_A1"],
#   clustervar = "walktrap_PCA_50",
#   return_plots = TRUE
#   # pdf_file = here("plots/01_build_spe/test_walk_trap_50.pdf")
# )
# 
# vis_grid_clus(spe,
#               clustervar = "walktrap_PCA_10",
#               pdf_file = here("plots/01_build_spe/test_walk_trap_10.pdf")
# )
# 


# plotReducedDim(spe, dimred="PCA", ncomponents=4,
#                color_by = "label")
# 
# plotReducedDim(spe, dimred="PCA", ncomponents=4,
#                color_by = "slide_num")


## TSNE --------------------------------------------------------------------
set.seed(100)
spe_all <- runTSNE(spe_all, dimred="PCA", perplexity=20)
# saveRDS(spe_all, "~/spe_pca.rds")




## UMAP --------------------------------------------------------------------
spe_all <- runUMAP(spe_all, dimred = "PCA") # TODO: edit this
# Same as below
# # set.seed(1)
# spe <- runUMAP(spe, dimred = "PCA", n_dimred = 10, name = "UMAP_10PCs_auto")
# spe_all <- runUMAP(spe_all, dimred = "PCA_10", n_dimred = 10, name = "UMAP_10PCs")# TODO: edit this
# spe_all <- runUMAP(spe_all, dimred = "PCA_10", n_dimred = 5, name = "UMAP_5PCs")# TODO: edit this

# ggpubr::ggarrange(
#   plotReducedDim(spe, dimred="UMAP_10PCs_auto", color_by = "sample_id"),
#   plotReducedDim(spe, dimred="UMAP_10PCs", color_by = "sample_id")
# )



# ggpubr::ggarrange(
#   plotUMAP(spe, color_by = "sample_id"), 
#   plotReducedDim(spe, dimred="UMAP_10PCs", color_by = "sample_id"),
#   common.legend = TRUE,
#   legend = "bottom"
# )
# 
# ggpubr::ggarrange(
#   plotUMAP(spe, color_by = "sample_id"), 
#   plotUMAP(spe, color_by = "walktrap_PCA_50"),
#   # common.legend = TRUE,
#   legend = "bottom"
# )
# 
# plotReducedDim(spe, dimred="UMAP_5PCs", color_by = "sample_id")
# 
# 
# plotUMAP(spe[, spe$sample_id==c("V12F14-057_A1", "V12F14-057_C1")],
#          color_by = "sample_id")
# plotUMAP(spe, color_by = "slide_num")
# plotUMAP(spe, color_by = "dx")


# Curated Genes  -------------------------------------------------------
## Select spatialDLPFC marker genes  ----------------------------------
file_DLPFC_enrich_csv <- here(
  "code/spatial_clustering/PRECAST",
  "TableS8_sig_genes_FDR5perc_enrichment.csv"
)

gene_df_raw <- read.csv(
  file_DLPFC_enrich_csv
)

n_marker_gene <- 100

gene_df <- gene_df_raw |> 
  filter(spatial_domain_resolution == "Sp09") |> 
  group_by(test) |> 
  arrange(fdr, .by_group = TRUE) |> 
  slice_head(n=n_marker_gene)


stopifnot(all(gene_df$model_type == "enrichment"))
stopifnot(nrow(gene_df) == 9*n_marker_gene)

cat("NOTE (boyiguo1): ",
    gene_df$ensembl |> unique() |> length(),
    " unique gene markers are selected for spatial clustering. \n")

spe_curated <- spe[gene_df$ensembl |> unique(), ]
spe_curated <- runPCA(spe_curated)

# TODO: create a variable for weird sample

spe_curated$weird_sample <- factor(spe_curated$sample_id) |> 
  forcats::fct_collapse(
    other_levels = setdiff(
      unique(spe_curated$sample_id), 
      c("V13F27-294_B1", "V12F14-057_A1")
    )
  )

# saveRDS(spe_curated, "~/spe_curated_pca.rds")
# spe_curated <- readRDS("~/spe_curated_pca.rds")
plot1 <- plotReducedDim(
  spe_curated, dimred="PCA", ncomponents=4,
  color_by = "sample_id"
) +
  theme(legend.position = "none")

plot1 |> ggsave(here("plots/spatial_cluster/test_pca_curated.pdf"))

# metadata(spe_curated)$dx_df
spe_curated$dx <- metadata(spe_curated)$dx_df$dx[
  match(
    spe_curated$sample_id,
    metadata(spe_curated)$dx_df$sample_id
  )]

# TODO: replace with esheR to make the point size much smaller
# TODO: play with gg merge R
plotReducedDim(
  spe_curated, dimred="PCA", ncomponents=2,
  color_by = "dx"
)


plotReducedDim(
  spe_curated, dimred="PCA", ncomponents=2,
  color_by = "weird_sample"
)

## TSNE --------------------------------------------------------------------
# Takes really long time to run
# spe_curated <- runTSNE(spe_curated, dimred="PCA", perplexity=20)
# saveRDS(spe_all, "~/spe_pca.rds")




## UMAP --------------------------------------------------------------------
spe_curated <- runUMAP(spe_curated, dimred = "PCA")
plotReducedDim(
  spe_curated, dimred="UMAP",
  color_by = "sample_id") +
  theme(legend.position = "none")

plotReducedDim(
  spe_curated, dimred="UMAP",
  color_by = "dx") +
  theme(legend.position = "none")
# TODO: edit this

# Session Info ------------------------------------------------------------
sessioninfo::session_info()

