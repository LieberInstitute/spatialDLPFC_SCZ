# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(pheatmap)
  library(sessioninfo)
})

# Read Data ----
## Dx DEG genes ----
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
)

# adj_p_cutoff <- 0.05
# sig_gene_df <- gene_df |>
#   filter(fdr_ntc <= adj_p_cutoff)

sig_gene <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    "Test_90DEGs.xlsx"
  ),
  col_names = FALSE
)[[1]]

sig_gene_df <- gene_df |>
  filter(gene %in% sig_gene)

# Annotation data for pheatmap
ann_df <- sig_gene_df |>
  column_to_rownames(var = "gene") |>
  transmute(
    SCZ_reg = factor(
      logFC_scz > 0,
      levels = c(TRUE, FALSE),
      labels = c("Up", "Down")
    )
  )

## PEC SpD Enrichment ----
# Local path
load(
  here(
    "code/analysis/pseudobulk_dx",
    "spatialDLPFC_BayesSpace_k09_registration.Rdata"
  )
)

pec_spd_df <- modeling_results$enrichment # |

spd_name_df <- read_csv(
  here("code/analysis/visium_spatial_clustering/bayesSpace_layer_annotations.csv")
)
# filter(spatial_domain_resolution == "Sp09") #


# pec_spd_gene_sub <- pec_spd_df |>
#   filter(pec_spd_df$ensembl %in% sig_gene_df$ensembl)


# Heatmap
# Orientation: gene (row) by spd (columns) with fill (-log10(p), maybe a directionality of logFC)

rownames(pec_spd_df) <- NULL

heatmap_pec_spd_df <- pec_spd_df |>
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

colnames(heatmap_pec_spd_df) <- spd_name_df$layer_combo[match(
  colnames(heatmap_pec_spd_df),
  spd_name_df$cluster
)]
# rownames(heatmap_pec_spd_df)

# stopifnot()

# mutate(
#   neg_log10_p = -1 * log10(pval)
# ) |>
# pivot_wider(
#   id_cols = c("ensembl"),
#   names_from = "test",
#   values_from = "neg_log10_p"
# ) |>
# column_to_rownames(var = "ensembl") |>
# select(starts_with("Sp09")) |>



# Make sure the order of genes is consistent with PEC spd.
gene_names_hc_ordered <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)


# saveRDS(
#   hc,
#   here("code/analysis/pseudobulk_dx",
#   "spd_hierarchical_cluster_order.rds")
# )

pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_pec_spd_p_value_centered.pdf"
  ),
  height = 20
)
heatmap_pec_spd_df[gene_names_hc_ordered, c(
  "Sp09D01 ~ L1",
  "Sp09D02 ~ L1",
  "Sp09D03 ~ L2",
  "Sp09D05 ~ L3",
  "Sp09D08 ~ L4",
  "Sp09D04 ~ L5",
  "Sp09D07 ~ L6",
  "Sp09D06 ~ WM",
  "Sp09D09 ~ WM"
)] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    # cluster_rows = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )
dev.off()




# Plot t-testistics ----

t_stat_mat <- pec_spd_df |>
  mutate(
    across(
      starts_with("t_stat"),
      ~ abs(.x),
      .names = "abs_{.col}"
    )
  ) |>
  filter(ensembl %in% sig_gene_df$ensembl) |>
  column_to_rownames(var = "gene") |>
  select(starts_with("abs_t_stat_")) |>
  rename_with(
    .fn = ~ gsub("^abs_t_stat_", "", .),
    .cols = starts_with("abs_t_")
  )

colnames(t_stat_mat) <- spd_name_df$layer_combo[
  match(colnames(t_stat_mat), spd_name_df$cluster)
]


pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_pec_spd_t_stat.pdf"
  ),
  height = 20
)
t_stat_mat[gene_names_hc_ordered, c(
  "Sp09D01 ~ L1",
  "Sp09D02 ~ L1",
  "Sp09D03 ~ L2",
  "Sp09D05 ~ L3",
  "Sp09D08 ~ L4",
  "Sp09D04 ~ L5",
  "Sp09D07 ~ L6",
  "Sp09D06 ~ WM",
  "Sp09D09 ~ WM"
)] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    # scale = "row",
    # cluster_rows = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )
dev.off()

pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_pec_spd_t_stat_row_centered.pdf"
  ),
  height = 20
)
t_stat_mat[gene_names_hc_ordered, c(
  "Sp09D01 ~ L1",
  "Sp09D02 ~ L1",
  "Sp09D03 ~ L2",
  "Sp09D05 ~ L3",
  "Sp09D08 ~ L4",
  "Sp09D04 ~ L5",
  "Sp09D07 ~ L6",
  "Sp09D06 ~ WM",
  "Sp09D09 ~ WM"
)] |>
  data.matrix() |>
  pheatmap(
    mat = _,
    scale = "row",
    # cluster_rows = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    annotation_row = ann_df
  )

dev.off()



## PEC snRNA Enrichment ----

# /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
# here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_velm_registration.Rdata")


load(
  # TODO: this is the wrong dataset, needs update after confirming with Louise
  # here(
  #   "sn_velm_registration.Rdata"
  # )
  here(
    "code/analysis/pseudobulk_dx/",
    "spatialDLPFC_sn_hc_registration.Rdata"
  )
)

pec_snRNA_df <- sn_hc_registration$enrichment

rownames(pec_snRNA_df) <- NULL

pdf(
  file = here(
    "plots/PB_dx_genes/",
    "test_sig_gene_enrich_pec_snRNA_p_value.pdf"
  ),
  height = 20
)
pec_snRNA_df |>
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
  ) |>
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
