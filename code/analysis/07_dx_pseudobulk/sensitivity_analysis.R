# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SpatialExperiment)
  # library(scater)
  # library(ggrepel)
  # library(spatialLIBD)
  # library(ggplot2)
  # library(ggrepel)
  library(sessioninfo)
})

# Load data ----
## Load data from prelim analysis ----
sce_prelim <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

## Load new data -----
sce_pseudo <- readRDS(
  here(
    "processed-data/rds/07_dx_pseudobulk",
    "sce_pseudo_PRECAST07_donor_spd.rds"
  )
)


# Examine difference in dimension ----
dim(sce_prelim) # old data
dim(sce_pseudo) # new data

# NOTE: The new data has less cells as well as less genes than the old data
# I don't understand why there's a difference in these numbers


# Old data has more of genes
nrow(sce_prelim) - nrow(sce_pseudo)
# [1] 3010

## These genes are
diff_gene_emsbl <- setdiff(rownames(sce_prelim), rownames(sce_pseudo))
rowData(sce_prelim)[diff_gene_emsbl, ]

gene_172_df <- read_csv(
  here("code/analysis/10_dx_deg_adjust_spd/172_prelim_fdr010.csv")
)

# Which previous 172 genes are missing in the new dat
# NOTE: there are some important genes are not included in the new pseudobulked data, such as "BDNF", "TEME119"
intersect(gene_172_df$ensembl, diff_gene_emsbl)
gene_172_df[gene_172_df$ensembl %in% diff_gene_emsbl, ] |>
  pull(gene)


## What are overletting genes





# Old data ocntains
# Old data has more of genes
ncol(sce_prelim) - ncol(sce_pseudo)

setdiff(
  paste0(sce_prelim$sample_id, ".", sce_prelim$PRECAST_07),
  paste0(sce_pseudo$sample_id, ".", sce_pseudo$PRECAST_07)
)



# Check difference in ncellss between the two datasets ----

tmp_data <- full_join(
  sce_prelim |> colData() |> data.frame() |> select(sample_id, PRECAST_07, ncells),
  sce_pseudo |> colData() |> data.frame() |> select(sample_id, PRECAST_07, ncells),
  by = c("sample_id", "PRECAST_07"),
  suffix = c("_prelim", "_pseudo")
)

tmp_data |> filter(ncells_prelim != ncells_pseudo)

# identical(
#   counts(sce_prelim)["ENSG00000157005", ],
#   counts(sce_pseudo)["ENSG00000157005", ]
# )

# Session Info ----
sessioninfo::session_info()
