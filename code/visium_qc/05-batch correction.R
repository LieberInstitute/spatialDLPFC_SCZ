# Require more than 10 GB, running with 80GB

# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(scran)
library(scater)
library(BiocParallel)
library(PCAtools)
library(harmony)


path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

spe <- readRDS(
  path_spe_after_spot_qc
)




# Create logcounts
spe <- logNormCounts(spe)



## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
                    block = spe$sample_id#,
                    # BPPARAM = MulticoreParam(4)
)

pdf(
  # TODO: change path
  here::here("plots", "01_build_spe", "test_scran_modelGeneVar_final.pdf"),
  useDingbats = FALSE
)
mapply(function(block, blockname) {
  plot(
    block$mean,
    block$total,
    xlab = "Mean log-expression",
    ylab = "Variance",
    main = blockname
  )
  curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
  )
}, dec$per.block, names(dec$per.block))
dev.off()

top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)
# [1] 2059

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 7478 instead of 18417 for DLPFC

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 5834 instead of 17830 for spatialDLPFC

# save(top.hvgs,
#      top.hvgs.fdr5,
#      top.hvgs.fdr1,
#      # TODO: change path
#      file = here::here("processed-data", "rdata", "spe", "01_build_spe", "top.hvgs_all.Rdata")
# )

set.seed(020122)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs.fdr5 , ncomponents = 50)
Sys.time()

dim(reducedDim(spe, "PCA"))


# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
# Currently 5 which is definitely wrong
# Compared to DLPFC
# 50

pdf(
  # TODO: edit path
  here::here("plots", "01_build_spe", "test_pca_elbow_final.pdf"),
  useDingbats = FALSE
)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()


summary(apply(reducedDim(spe, "PCA"), 2, sd))

# RunUMAP
set.seed(030122)
Sys.time()
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()

# Run TSNE
# set.seed(030122)
# Sys.time()
# spe <-
#   runTSNE(spe,
#           dimred = "PCA",
#           name = "TSNE_perplexity80",
#           perplexity = 80
#   )
# Sys.time()


# make plots of UMAP
# TODO: replace code with ggblend
# pdf(file = here::here("plots", "01_build_spe", "UMAP_subject.pdf"))
# ggplot(
#   data.frame(reducedDim(spe, "UMAP")),
#   aes(x = UMAP1, y = UMAP2, color = factor(spe$subject))
# ) +
#   geom_point() +
#   labs(color = "Subject") +
#   theme_bw()
# dev.off()

pdf(file = here::here("plots", "01_build_spe", "test_UMAP_sample_id.pdf"))
ggplot(
  data.frame(reducedDim(spe, "UMAP")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()


### harmony batch correction
spe <- RunHarmony(spe, "sample_id", verbose = F)
set.seed(1)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

# TODO: change path
pdf(file = here::here("plots", "01_build_spe", "test_UMAP_harmony_sample_id.pdf"))
# Oh, my god, perfect match?
# TODO: check the presentation for a subset of samples that shares
#      different composition of morphology.
ggplot(
  data.frame(reducedDim(spe, "UMAP.HARMONY")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()


########### for bayesSpace

## do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
# TODO: rename the row and col to show it is concatenated.
spe$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- colData(spe)$array_col

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46

# TODO: to remove this after generating new data. doesn't needs this because dx_df becomes a element
metadata(spe) <- list()
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")

# Sys.time()
#TODO: change path
save(spe, file = here::here("processed-data", "rds", "spe",
                            "01_build_spe", "test_spe_filtered_final.rds"))
# Sys.time()