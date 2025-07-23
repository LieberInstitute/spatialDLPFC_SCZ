# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  library(limma)
  library(sessioninfo)
  library(here)
})

# Load data ----
## Load pseudobulked data ----
sce_pseudo <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_PRECAST07_donor_spd.rds"
  )
)

## Load PRS dat ----
prs_data <- read_csv(
  here(
    "processed-data/donor_prs",
    "Spatial_DLPFC_SCZ_PRS.csv"
  )
) |>
  rename_with(
    ~ paste0("genotype_", .x), starts_with("PC")
  ) |>
  select(
    brnum = IID,
    PRS,
    starts_with("genotype_")
  )


colData(sce_pseudo) <- as.data.frame(colData(sce_pseudo)) |>
  left_join(prs_data, by = "brnum") |>
  DataFrame()

# limma test ----
## Create desgin matrix ----
dx_mod <- model.matrix(
  ~ 0 + PRS + fnl_spd + age + sex + slide_id + genotype_PC1 + genotype_PC2,
  colData(sce_pseudo)
)

## Calculate correlation within block ----
corfit <- duplicateCorrelation(
  logcounts(sce_pseudo),
  design = dx_mod,
  block = colData(sce_pseudo)$sample_id
)
## Run limma fit ----
fit <- lmFit(
  logcounts(sce_pseudo),
  design = dx_mod,
  block = colData(sce_pseudo)$sample_id,
  correlation = corfit$consensus
)

## Fit empirical Bayes model ----
eb_fit <- eBayes(fit)

test_res <- topTable(eb_fit,
  coef = "PRS", number = Inf, genelist = fit$genes,
  adjust.method = "BH", sort.by = "B", resort.by = NULL,
  p.value = 1, fc = NULL, lfc = NULL, confint = FALSE
)

fnl_res <- test_res |>
  rownames_to_column("gene_id") |>
  left_join(
    rowData(sce_pseudo) |>
      as.data.frame() |>
      select(gene_id, gene_name)
  )

# Save results ----
write_csv(fnl_res, file = here(
  "processed-data/rds/14_prs_deg",
  "PRS_DEG_test_res_PRECAST07_donor_spd.csv"
))


## TEST (deprecated): Use SpatialLIBD wrapper ----
## NOTE: SpatialLIBD wrapper would give inaccurate results due to technical reasons
## althought using the SpatialLIBD Wrapper can give you a result.
## Technically, the result is inaccurate, eecasue PRS is not categorical, and hence the wrapper function doesn't function with continuous varibale properly.
# dx_mod <- registration_model(
#   sce_pseudo,
#   covars = c("fnl_spd", "age", "sex", "slide_id", "genotype_PC1", "genotype_PC2"),
#   var_registration = "PRS"
# )

# dx_block_cor <- registration_block_cor(
#   sce_pseudo,
#   registration_model = dx_mod,
#   var_sample_id = "sample_id"
# )

# dx_res <- registration_stats_enrichment(
#   sce_pseudo,
#   block_cor = dx_block_cor,
#   covars = c("fnl_spd", "age", "sex", "slide_id", "genotype_PC1", "genotype_PC2"),
#   var_registration = "PRS",
#   gene_ensembl = "gene_id",
#   gene_name = "gene_name"
# )

# Session Info ----
session_info()
