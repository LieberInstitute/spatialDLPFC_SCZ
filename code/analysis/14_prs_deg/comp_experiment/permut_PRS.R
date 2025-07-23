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

# Create design matrix ----

colData(sce_pseudo) <- as.data.frame(colData(sce_pseudo)) |>
  left_join(prs_data, by = "brnum") |>
  DataFrame()



# Permute PRS column ----
set.seed(123)  # For reproducibility

sce_pseudo_permuted <- sce_pseudo
sce_pseudo_permuted$PRS <- sample(sce_pseudo_permuted$PRS)


# Fit model ----
## Create desgin matrix ----
dx_mod <- model.matrix(
  ~ 0 + PRS + fnl_spd + age + sex + slide_id + genotype_PC1 + genotype_PC2,
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
eb_fit <- eBayes(fit)

test_res <- topTable(eb_fit,
  coef = "PRS", number = Inf, genelist = fit$genes,
  adjust.method = "BH", sort.by = "B", resort.by = NULL,
  p.value = 1, fc = NULL, lfc = NULL, confint = FALSE
)

test_res$P.Value |> hist()

volcanoplot(eb_fit,
  coef = 1, style = "p-value",
  highlight = 10, names = fit$genes$ID,
  hl.col = "blue", xlab = "Log2FC (PRS)", ylab = NULL, pch = 16, cex = 0.35
)

ggplot() +
  aes(
    x = test_res$logFC,
    y = -log10(test_res$adj.P.Val)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_point(alpha = 0.5) +
  labs(
    x = "Log2 Fold Change (PRS)",
    y = "-log10(FDR adjusted P-value)",
    title = "Volcano Plot: PRS DEGs"
  ) +
  theme_minimal()


adj_sig_DEG_name <- test_res |> filter(
  adj.P.Val < 0.05,
) |> rownames()

rowData(sce_pseudo)[adj_sig_DEG_name, "gene_name"]

# Save results ----


# Session Info ----
session_info()
