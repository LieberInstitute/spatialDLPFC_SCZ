# JHPCE Interactive -------------
# srun --pty --x11 --mem=80G bash


# Load Packages --------------------
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(scater)
  library(scran)
  # library(pryr)                 # Check spe size
  # library(spatialLIBD)
  library(tidyverse)
})


# Load raw data ------
path_raw_spe <- here(
  "processed-data/rds",
  "01_build_spe",
  "raw_spe_wo_SPG_N63.rds"
)

raw_spe <- readRDS(
  path_raw_spe
)

# Remove spots with 0 transcripts
raw_spe <- raw_spe[, raw_spe$sum_umi != 0]

# Create logcounts
if (!"logcounts" %in% assayNames(raw_spe)) {
  raw_spe <- logNormCounts(raw_spe)
}



# Selected Genes ------------------------------------
## Highly Variable Genes ----
raw_spe <- scran::modelGeneVar(sce)

# TODO: what is the different between block and non-block?

# TODO: is there a way to subseting the genes in runPCA

## spatialDLPFC marker genes ----------------------------------
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
  slice_head(n = n_marker_gene)

# Error Prevention
stopifnot(all(gene_df$model_type == "enrichment"))
stopifnot(nrow(gene_df) == 9 * n_marker_gene)

cat(
  "NOTE (boyiguo1): ",
  gene_df$ensembl |> unique() |> length(),
  " unique gene markers are selected for spatial clustering. \n"
)

genes_spatialDLPFC <- gene_df$ensembl |> unique()





# Dim Red & Viz ----
## PCA ----
spe_curated <- runPCA(spe[genes_spatialDLPFC, ])

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
  spe_curated,
  dimred = "PCA", ncomponents = 4,
  color_by = "sample_id"
) +
  theme(legend.position = "none")

plot1 |> ggsave(here("plots/spatial_cluster/test_pca_curated.pdf"))

# metadata(spe_curated)$dx_df
spe_curated$dx <- metadata(spe_curated)$dx_df$dx[
  match(
    spe_curated$sample_id,
    metadata(spe_curated)$dx_df$sample_id
  )
]

# TODO: replace with esheR to make the point size much smaller
# TODO: play with gg merge R
plotReducedDim(
  spe_curated,
  dimred = "PCA", ncomponents = 2,
  color_by = "dx"
)


plotReducedDim(
  spe_curated,
  dimred = "PCA", ncomponents = 2,
  color_by = "weird_sample"
)

## TSNE --------------------------------------------------------------------
# Takes really long time to run
# spe_curated <- runTSNE(spe_curated, dimred="PCA", perplexity=20)
# saveRDS(spe_all, "~/spe_pca.rds")




## UMAP --------------------------------------------------------------------
spe_curated <- runUMAP(spe_curated, dimred = "PCA")
plotReducedDim(
  spe_curated,
  dimred = "UMAP",
  color_by = "sample_id"
) +
  theme(legend.position = "none")

plotReducedDim(
  spe_curated,
  dimred = "UMAP",
  color_by = "dx"
) +
  theme(legend.position = "none")
# TODO: edit this

# Session Info -----------------------
sessioninfo::session_info()
