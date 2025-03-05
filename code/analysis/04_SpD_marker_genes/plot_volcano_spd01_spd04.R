# Load library ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
  library(sessioninfo)
})

# Load data ----
## Read pairwise test res ----
layer_res <- readRDS(here(
  # TODO: need organize
  "processed-data", "rds", "layer_enrich_test",
  "test_pairwise_PRECAST_07.rds"
)) |>
  select(ends_with("spd01-spd04"), ensembl, gene) |>
  rename_with(~ gsub("_spd01-spd04", "", .))


## Read Sang Ho's annotation ----
df <- read_csv(
  here(
    "code/analysis/04_SpD_marker_genes",
    "pairwise_DEG_enriched_in_spd01_highlights.csv"
  ),
  col_names = FALSE
) |>
  select(
    "type" = X1,
    "gene_name" = X2
  ) |>
  mutate(
    gene_name = str_trim(gene_name)
  )

layer_res <- layer_res |> left_join(
  df,
  by = join_by(gene == gene_name)
)

table(layer_res$type)

stopifnot(sum(table(layer_res$type)) == nrow(df))
# layer_res |> filter(type == "Glutamatergic")
# layer_res |> filter(gene_name == "GRIN1")


# Make volcano plot ----
volcano_p <- ggplot() +
  geom_point(
    data = layer_res |> filter(is.na(type)),
    aes(x = logFC, y = -log10(fdr), color = type)
  ) +
  geom_point(
    data = layer_res |> filter(!is.na(type)),
    aes(x = logFC, y = -log10(fdr), color = type)
  ) +
  geom_hline(aes(yintercept = -log10(0.05))) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("FDR")),
    color = "Type"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Save plot ----
ggsave(
  here(
    "plots/04_SpD_marker_genes",
    "volcano_plot_spd01-04.pdf"
  ),
  volcano_p,
)

# Session info ----
sessioninfo::session_info()
