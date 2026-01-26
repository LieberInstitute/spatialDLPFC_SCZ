# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(limma)
  library(scater)
  library(tidyverse)
  library(sessioninfo)
})

# Load pseudobulk data ----
sce_pseudo <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_PRECAST07_donor_spd.rds"
  )
)

# Interaction model desgin matrix ----
dx_mod <- model.matrix(
  ~ 0 + dx * fnl_spd + age + sex + slide_id,
  colData(sce_pseudo)
)


# Run limma fit ----
## Calculate correlation within block ----
corfit <- duplicateCorrelation(
  logcounts(sce_pseudo),
  design = dx_mod,
  block = colData(sce_pseudo)$sample_id
)

fit <- lmFit(
  logcounts(sce_pseudo),
  design = dx_mod,
  block = colData(sce_pseudo)$sample_id,
  correlation = corfit$consensus
)

## Fit empirical Bayes model ----
fit <- eBayes(fit)

# Save results ----
saveRDS(
  fit,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "limma_obj_int_PRECAST_07_redo.rds"
  )
)

# Session Info ----
sessioninfo::session_info()
