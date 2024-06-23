# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(pheatmap)
  library(sessioninfo)
})

# Read Data ----
.dx_group <- "ntc"
# .dx_group <- "scz"

## Dx DEG genes ----
gene_df <- read_csv(
    here(
      "processed-data/PB_dx_genes/",
      "test_PRECAST_07.csv"
    )
  # )

  # "~/Downloads/test_PRECAST_07.csv"
)

sig_gene_df <- gene_df |>
  filter(fdr_ntc <= 0.05)


ann_df <- gene_df |>
  filter(fdr_ntc <= 0.05) |>
  column_to_rownames(var = "gene") |>
  transmute(
    SCZ_reg = factor(
      logFC_scz > 0,
      levels = c(TRUE, FALSE),
      labels = c("Up", "Down")
    )
  )

spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

# Enrichment test results ----
modeling_results <- readRDS(
  here(
    "processed-data", "rds", "layer_enrich_test",
    paste0("test_enrich_", .spd, "_", .dx_group, ".rds")
  )
)

rownames(modeling_results) <- NULL

heatmap_pec_spd_df <- modeling_results |>
  mutate(
    across(
      starts_with("p_value"),
      ~ -1 * log10(.x),
      .names = "neg_log10_{.col}"
    )
  ) |>
  filter(ensembl %in% sig_gene_df$ensembl) |>
  column_to_rownames(var = "gene") |>
  select(starts_with("neg_log10_")) |>
  rename_with(
    .fn = ~ gsub("^neg_log10_p_value_", "", .),
    .cols = starts_with("neg_log10_")
  )

colnames(heatmap_pec_spd_df) <- spd_anno_df$anno_lab[
  match(colnames(heatmap_pec_spd_df), spd_anno_df$spd)
]



pdf(
  file = here(
    "plots/PB_dx_genes/",
    paste0("test_sig_gene_enrich_SCZ_spd_",.dx_group,".pdf")
  ),
  height = 20
)
heatmap_pec_spd_df[, order(colnames(heatmap_pec_spd_df))] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )
dev.off()

# Session Info ----
session_info()
