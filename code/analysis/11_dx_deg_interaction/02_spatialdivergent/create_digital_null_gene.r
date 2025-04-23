suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
})


# Create Permuted Pseudobulked data ----
set.seed(20240724)

## Load pseudobulk data ----
sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

sce_pseudo$spd <- sce_pseudo$PRECAST_07

## subsample the dataset ----
all_genes <- rownames(sce_pseudo)
rand_genes <- sample(all_genes, size = 100) # 100 random genes

# Subset genes
tmp_pseudo <- sce_pseudo[rand_genes, ]



## Permute in logcount_matrix data matrix

# log_mat <- logcounts(tmp_pseudo)

# Create permuted column order - perm spd within each individual
# NOTE: optimally, I should permute the spd within each indiviudal for each gene. however, it is technically challenge to do.

col_df <- colData(tmp_pseudo) |> data.frame()

perm_spd_df <- col_df |>
  # group_by(sample_id) |>
  reframe(
    # n = n(),
    spd = spd,
    perm_spd = sample(spd),
    .by = "sample_id"
  )

# create random order per sample
# perm_spd_df <- map2_dfr(
#   n_spd_per_sample$sample_id,
#   n_spd_per_sample$n,
#   .f = function(.sample_id, .n) {
#     data.frame(
#       sample_id = .sample_id,
#       raw_spd = sprintf("spd%02d", seq.int(1, .n)),
#       perm_spd = sample( # Permute the label
#         sprintf("spd%02d", seq.int(1, .n))
#       )
#     )
#   }
# )

# Left_join to the col_df
col_df <- col_df |>
  left_join(
    perm_spd_df,
    by = join_by(
      sample_id == sample_id,
      spd == spd
    ),
    unmatched = "error"
  )

colData(tmp_pseudo) <- col_df |> DataFrame()


# Fit the interaction model ----
# TODO: copy form create_limma_fit_interaction.R and hypo_test_dx_deg_across_spd.r
dx_mod <- model.matrix(
  ~ 0 + dx * perm_spd + age + sex + slide_id,
  colData(tmp_pseudo)
)

int_terms <- grep(":", dx_mod |> colnames(), value = TRUE)

# Run limma fit ----
## Calculate correlation within block
corfit <- duplicateCorrelation(
  logcounts(tmp_pseudo),
  design = dx_mod,
  block = colData(tmp_pseudo)$sample_id
)

fit <- lmFit(
  logcounts(tmp_pseudo),
  design = dx_mod,
  block = colData(tmp_pseudo)$sample_id,
  correlation = corfit$consensus
)
fit <- eBayes(fit)

cont.mat <- rbind(
  rep(-1, 7),
  rep(1, 7),
  matrix(0, nrow = 23, ncol = 7),
  cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
)
colnames(cont.mat) <- sprintf("spd%02d", 1:7)

# Test any non-zero FC across SpD ----
contrast_fit <- contrasts.fit(fit, cont.mat)
contrast_fit <- eBayes(contrast_fit)

cont_df <- topTable(contrast_fit,
  coef = sprintf("spd%02d", 1:7),
  num = Inf
) |>
  rownames_to_column(var = "ensembl") |>
  mutate(
    ensembl = paste0("perm_", ensembl),
    gene = sprintf("perm_gene%03d", seq.int(1, n()))
  )


# Validation against non-permuted data
# log2fc_mat <- readRDS(
#   here(
#     "processed-data/PB_dx_genes/interaction",
#     "test_laminar_specific_log2FC.rds"
#   )
# )

# log2fc_mat |> filter(ensembl == "ENSG00000064726")

## Save logFC estimates from 
cont_df |>
  select(
    ensembl, gene,
     starts_with("spd")
    ) |> 
saveRDS(
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_permed_laminar_specific_log2FC.rds"
  )
)



# Session Info ----
sessioninfo::session_info()
