library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)
library(scran)
library(scater)
library(BiocParallel)
library(PCAtools)

path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

spe <- readRDS(
  path_spe_after_spot_qc
)


# Create logcounts
stopifnot("logcounts" %in% assayNames(spe))


## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
                    block = spe$sample_id#,
                    # BPPARAM = MulticoreParam(4)
)

# Visualization of TODO: ? ---------
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

fld_ftr_slct <- here("processed-data/rds/feature_selection")

dir.create(fld_ftr_slct,
           recursive = TRUE, showWarnings = FALSE)

saveRDS(dec, file = here(fld_ftr_slct, "HVG_imp2.rds"))

# HVG subsetting ======
top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)
# [1] 2059

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 7478 instead of 18417 for DLPFC

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)