suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
})


# Create Permuted Pseudobulked data ----
set.seed(20240719)

## Load pseudobulk data ----
sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

## subsample the dataset ----
all_genes <- rownames(sce_pseudo)
rand_genes <- sample(all_genes, size = 100) # 100 random genes

# Subset genes
tmp_pseudo <- sce_pseudo[rand_genes]

## Permute in logcount_matrix data matrix

log_mat <- logcounts(tmp_pseudo)

# Create permuted column order - perm spd within each individual
# NOTE: optimally, I should permute the spd within each indiviudal. however, it is technically challenge to do.

col_df <- colData(sce_pseudo) |> data.frame()

n_spd_per_sample <- col_df |>
  group_by(sample_id) |>
  summarize(n = n())

# create random order per sample
perm_spd_df <- n_spd_per_sample |> map2_dfr(
  .f = function(.sample_id, .n) {
    data.frame(
      sample_id = .sample_id,
      raw_spd = sprintf("spd_%02d", seq.int(1, .n)),
      perm_spd = sample( # Permute the label
        sprintf("spd_%02d", seq.int(1, .n))
      )
    )
  }
)

# Left_join to the col_df
col_df <- col_df |>
  left_join(
    perm_spd_df,
    by = c(
      "sample_id" = "sample_id",
      "spd" = "raw_spd"
    )
  )


# Fit the interaction model ----
# TODO: copy form create_limma_fit_interaction.R and hypo_test_dx_deg_across_spd.r





# Session Info ----
sessioninfo::session_info()
