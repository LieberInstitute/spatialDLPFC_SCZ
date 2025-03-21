# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(limma)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(pheatmap)
})

# Configuration ----
.file <- "test_spe_pseudo_PRECAST_07.rds"

.spd <- str_remove(.file, "test_spe_pseudo_") |>
  str_remove(".rds")


# Load pseudobulk data ----
sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    .file
  )
)

sce_pseudo$spd <- sce_pseudo[[.spd]]

# Interaction model desgin matrix ----
dx_mod <- model.matrix(
  ~ 0 + dx * spd + age + sex + slide_id,
  colData(sce_pseudo)
)

int_terms <- grep(":", dx_mod |> colnames(), value = TRUE)

# Run limma fit ----
## Calculate correlation within block
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
fit <- eBayes(fit)

# Save results ----
saveRDS(
  fit,
  here(
    "processed-data/PB_dx_genes",
    "test_inter_PRECAST_07_20240627.rds"
  )
)

# Session Info ----
sessioninfo::session_info()
