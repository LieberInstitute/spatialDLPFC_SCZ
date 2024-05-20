# Load packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(sessioninfo)
})


# Load Data ----

spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

gene_df <- read_csv(
  #   here(
  #     "processed-data/PB_dx_genes/",
  #     "test_PRECAST_07.csv"
  #   )
  # )

  "~/Downloads/test_PRECAST_07.csv"
)

sig_gene_df <- gene_df |>
  filter(fdr_ntc <= 0.05)

# sig_genes_df <- sig_gene_df <- gene_df |>
#   filter(fdr_ntc <= 0.05)

# up_genes <- sig_genes_df

# down_genes <-




# Heat map ----

col_df <- colData(spe_pb) |>
  data.frame() |>
  select(PRECAST_07, dx, sample_id) |>
  mutate(
    PRECAST_07 = factor(PRECAST_07,
      levels = c(
        "spd07", "spd06", "spd02",
        "spd05", "spd03", "spd01", "spd04"
      )
    )
  ) |>
  arrange(
    PRECAST_07, dx, sample_id
  )

gene_mat <- logcounts(spe_pb)[sig_gene_df$ensembl, rownames(col_df)]


# Orignial scale ----
pheatmap(
  mat = gene_mat,
  scale = "none",
  cluster_cols = FALSE,
  annotation_col = col_df |> select(
    PRECAST_07,
    # sample_id, # Overwhelm color pallete
    dx
  ),
  # Turn off names
  show_rownames = FALSE,
  show_colnames = FALSE,
)

# Row Centered plot -----
centered_gene_mat <- apply(gene_mat, 1, scale, scale = FALSE) |> t()

colnames(centered_gene_mat) <- colnames(gene_mat)

pheatmap(
  mat = centered_gene_mat,
  scale = "none",
  cluster_cols = FALSE,
  annotation_col = col_df |> select(
    PRECAST_07,
    # sample_id, # Overwhelm color pallete
    dx
  ),
  # Turn off names
  show_rownames = FALSE,
  show_colnames = FALSE,
)

# Row centered and scaled ----
pheatmap(
  mat = gene_mat,
  scale = "row",
  cluster_cols = FALSE,
  annotation_col = col_df |> select(
    PRECAST_07,
    # sample_id, # Overwhelm color pallete
    dx
  ),
  # Turn off names
  show_rownames = FALSE,
  show_colnames = FALSE,
)

# Session Info ----
session_info()
