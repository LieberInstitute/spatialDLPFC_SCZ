# Load packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SpatialExperiment)
  library(Seurat)
  library(PRECAST)
  library(sessioninfo)
  library(readxl)
  library(ComplexHeatmap)
  library(viridis)
  library(escheR)
  library(aricode)
  library(Polychrome)
})


# Load Necessary Data ----
## Raw QC-ed data ----
raw_spe <- readRDS(
  here::here(
    "processed-data/rds/02_visium_qc",
    "qc_spe_wo_spg_N63.rds"
  )
)
### 20 Samples ----
# Read Excel file
sub_samples <- read_xlsx(
  here::here("code/xenium_panel_design/Xenium_DonorList_Edit.xlsx"),
  col_names = FALSE
)[, 1:2] |> unlist()
sub_samples <- sub_samples[!is.na(sub_samples)]
stopifnot(length(sub_samples) == 24)
spe <- raw_spe[, raw_spe$brnum %in% sub_samples]
stopifnot(length(unique(spe$sample_id)) == 24)


## Pseudobulked data
raw_spe_pb <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)
spe_pb <- raw_spe_pb[, raw_spe_pb$BrNumbr %in% sub_samples]



# Define Gene Panels -----
## TODO: add gene names heres
# raw_gene_df <- read_csv(
#   here(
#     "processed-data",
#     "test_xenium_design_layer_genes_deconvo_buddies_all.csv"
#   )
# ) |> mutate(
#   spd = cellType.target,
#   ensembl = gene_ensembl
# )

# gene_df <-  raw_gene_df |>
#   group_by(cellType.target) |>
#   slice_head(n = 10) |> ungroup()

# gene_names <- gene_df$ensembl


raw_gene_df <- read.csv(
  here(
    # TODO: change path
    "code/analysis/visium_spatial_clustering",
    "TableS8_sig_genes_FDR5perc_enrichment.csv"
  )
) |>
  mutate(spd = test)

gene_df <- raw_gene_df |>
  filter(spatial_domain_resolution == "Sp09") |>
  group_by(test) |>
  arrange(fdr, .by_group = TRUE) |>
  slice_head(n = 10)

gene_df <- gene_df[!gene_df$gene |> duplicated(), ]

gene_names <- gene_df$ensembl

gene_df |>
  ungroup() |>
  select(ensembl, gene) |>
  write_csv(here(
    "code/xenium_panel_design",
    "layer_gene_86.csv"
  ))

# Create Heatmap -----
source(
  here("code/xenium_panel_design/heatmap_layer_probes.r")
)

p_heatmap <- create_heatmap(spe_pb, gene_names)

# Run PRECAST -----
source(
  here("code/xenium_panel_design/layer_prob_via_PRECAST.R")
)

# NOTE: this step may take some time.
result_spe <- run_PRECAST_with_genes(spe, gene_names)


# Calculate concordance -----
# Qualitative Examination ----
# col_dat <- colData(result_spe) |> data.frame()

# col_dat |>
#   group_by(sample_id) |>
#   summarize(ari = ARI(PRECAST_07, PRECAST_07_test))

plot_concordance <- function(col_dat) {
  spd_anno_df <- read_csv(
    here(
      "processed-data/man_anno",
      "spd_labels_k7.csv"
    )
  ) |>
    mutate(anno_lab = paste0(label, " (", spd, ") "))

  col_dat |>
    group_by(sample_id, PRECAST_07) |>
    summarize(
      n = n(),
      concord = RI(PRECAST_07, PRECAST_07_test)
    ) |>
    left_join(
      metadata(spe)$dx_df |> select(sample_id, dx),
      by = "sample_id"
    ) |>
    ggplot() +
    # geom_point(aes(x = PRECAST_07, y = concord, color = dx, group = dx)) +
    geom_boxplot(aes(x = PRECAST_07, y = concord, color = dx)) +
    scale_x_discrete(
      limits = spd_anno_df$spd[order(spd_anno_df$anno_lab)],
      labels = setNames(spd_anno_df$anno_lab, spd_anno_df$spd)
    ) +
    theme_minimal()
}

# Create spot plot ----

# Compile plots to a file ----
pdf(
  here(
    "plots/xenium_panel_design",
    paste0("test_marker_gene_", length(gene_names), ".pdf")
    #  "_per_layer_deconvo_buddies.pdf")
  ),
  height = 4, width = 10
)
## Heatmap ----
p_heatmap |> print()

## Concordance plot ----
# TODO: customize file name
result_spe |>
  colData() |>
  data.frame() |>
  plot_concordance() |>
  print()


## Spot Plot ----
# TODO: customize file name

for (.sample_id in unique(result_spe$sample_id)) {
  ggpubr::ggarrange(
    make_escheR(result_spe[, result_spe$sample_id == .sample_id]) |>
      add_fill("PRECAST_07") +
      scale_fill_manual(
        values = set_names(
          Polychrome::palette36.colors(7)[seq.int(7)],
          unique(result_spe[["PRECAST_07"]]) |> sort()
        ),
        guide = guide_legend(title = NULL)
      ) + labs(title = .sample_id),
    make_escheR(result_spe[, result_spe$sample_id == .sample_id]) |>
      add_fill("PRECAST_07_test") +
      scale_fill_manual(
        values = set_names(
          Polychrome::palette36.colors(7)[seq.int(7)],
          unique(result_spe[["PRECAST_07"]]) |> sort()
        ),
        guide = guide_legend(title = NULL)
      ) + labs(title = .sample_id)
  ) |> print()
}
dev.off()

# Session info ----
session_info()
