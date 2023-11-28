# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)


path_spe_after_spot_qc <- here::here("processed-data", "rds",
                                     "spe", "spe_after_spot_qc.rds")

spe <- readRDS(
  path_spe_after_spot_qc
)



# Feature/Gene Analysis -----------------------------------------------------
# * Remove Genes with 0 count ---------------------------------------------

gene_sum <- rowSums(counts(spe))
sum(gene_sum == 0)
# [1] 7229
mean(gene_sum == 0)     # High proportion of non expressing gene
# [1] 0.1975083

spe <- spe[gene_sum != 0,]
stopifnot(all(rowSums(counts(spe))!=0))

# * (Optional) Remove Genes with 0 variance ---------------------------------------------
gene_var <- counts(spe) |> MatrixGenerics::rowVars()
if(mean(gene_var == 0) != 0){
  spe <- spe[gene_var != 0, ]
}

stopifnot(all(MatrixGenerics::rowVars(counts(spe))!=0))
