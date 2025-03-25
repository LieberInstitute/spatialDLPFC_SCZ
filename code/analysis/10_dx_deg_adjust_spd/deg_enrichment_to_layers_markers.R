# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(spatialLIBD)
  library(sessioninfo)
})

# Load Data ----
## Load DEGs ----
gene_df <- read_csv(
  here(
    "processed-data/rds/10_dx_deg_adjust_spd",
    "dx-deg_PRECAST07.csv"
  )
)

## format enrichment test res
t_stats <- gene_df[, grep("^t_stat_", colnames(gene_df))]
colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))

## Load layer-markers ----
raw_layer_df <- read_csv(
  here(
    "code/analysis/04_SpD_marker_genes",
    "raw_enrichment_spd_marker_gene.csv"
  )
)


## Load Nat Neruo layer-markers ----
load(
  here(
    "code/analysis/05_spatial_registration",
    "Human_DLPFC_Visium_modeling_results.Rdata"
  )
)

nat_neuro_layer_df <- modeling_results

## Load PEC layer-markers ----

## Format the layer_df to a list of layer-specific gene sets ----
layer_df_long <- raw_layer_df |>
  select(ensembl, gene, starts_with("fdr")) |>
  # turn to long form
  pivot_longer(
    cols = starts_with("fdr"),
    names_to = "layer",
    values_to = "fdr"
  ) |>
  mutate(
    layer = str_remove(layer, "fdr_")
  ) |>
  filter(
    fdr < 0.05
  )

list_layer_marker_genes <- layer_df_long |>
  group_by(
    layer
  ) |>
  summarise(
    ensembl = list(ensembl)
  ) |>
  deframe()



# Enrichment analysis ----
spatialLIBD::gene_set_enrichment(
  gene_list = list_layer_marker_genes,
  modeling_results = list("enrichment" = gene_df),
  model_type = "enrichment",
  fdr_cut = 0.05
)
# Unexpected observation: why the NumSig is different between ntc and scz tests?
# I would expect NumSig to be the same for ntc and scz tests.
# Is this a bug of the function?

# This is because that the only chhose positive logFC genes for scz test for each test
# I'm not sure if it is the best impelmentation for this type of analysis

# Check if fdr_scz and fdr_ntc columns are identical up to a reasonable precision level

identical_fdr <- all.equal(gene_df$fdr_scz, gene_df$fdr_ntc, tolerance = 1e-8)
print(identical_fdr)

gene_df[with(gene_df, t_stat_ntc == -1 * t_stat_scz), ]



res <- spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = list("enrichment" = raw_layer_df),
  model_type = "enrichment",
  fdr_cut = 0.05
)

#          OR         Pval  test NumSig SetSize       ID model_type fdr_cut
# 1 1.5740054 4.234266e-03 spd01     57     172 172 degs enrichment    0.05
# 2 0.6205032 9.983593e-01 spd02     47     172 172 degs enrichment    0.05
# 3 0.5134703 9.997538e-01 spd03     27     172 172 degs enrichment    0.05
# 4 1.6832864 5.186548e-04 spd04     81     172 172 degs enrichment    0.05
# 5 0.4941226 9.999815e-01 spd05     38     172 172 degs enrichment    0.05
# 6 0.5959785 9.990166e-01 spd06     42     172 172 degs enrichment    0.05
# 7 2.6533904 2.891710e-10 spd07     89     172 172 degs enrichment    0.05


# Make visualizaiton of the results ----
# Load annotation
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))

spd_order <- order(spd_anno_df$anno_lab)

res |>
  inner_join(
    spd_anno_df,
    by = c("test" = "spd")
  ) |>
  mutate(test = anno_lab) |>
  select(-label, -anno_lab) |>
  gene_set_enrichment_plot(
    # res,
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    cex = 1.5 # control the size of the text
  ) + title(
    "scz-DEGs enriched in PNN project"
  )

## Re-order and clean up the SpD annotaiton ----
## Enrich dx-DEGs in Nat Neuro ----
spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = nat_neuro_layer_df,
  model_type = "enrichment",
  fdr_cut = 0.05
) |> gene_set_enrichment_plot(
  PThresh = 12,
  ORcut = 3,
  enrichOnly = FALSE,
  cex = 1.5 # control the size of the text
) +
  title("SCZ-DEGs enriched in Nat Neruo")


## Enrich dx-DEG in PEC study

pec_layer_df <- fetch_data(
  type = "spatialDLPFC_Visium_modeling_results",
  destdir = tempdir(),
  eh = ExperimentHub::ExperimentHub(),
  bfc = BiocFileCache::BiocFileCache()
)

bayes_anno <-
  read_csv(
    file = file.path(
      "code/analysis/10_dx_deg_adjust_spd",
      "bayesSpace_layer_annotations.csv"
    )
  ) |>
  mutate(
    anno_lab = paste0(layer_annotation, " (", cluster, ") ")
  ) |>
  select(layer_combo,
    test = cluster,
    Annotation = anno_lab
  )

spatialLIBD::gene_set_enrichment(
  gene_list = list(
    `172 degs` = gene_df |> filter(fdr_scz < 0.10) |> pull(ensembl),
    `upreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz > 0) |> pull(ensembl),
    `downreg genes` = gene_df |> filter(fdr_scz < 0.10 & logFC_scz < 0) |> pull(ensembl)
  ),
  modeling_results = pec_layer_df,
  model_type = "enrichment",
  fdr_cut = 0.05
) |>
  left_join(
    bayes_anno,
    by = c("test" = "test")
  ) |>
  mutate(
    test = Annotation
  ) |>
  gene_set_enrichment_plot(
    PThresh = 12,
    ORcut = 3,
    enrichOnly = FALSE,
    cex = 1.5 # control the size of the text
  ) +
  title("SCZ-DEGs enriched in PEC")

# Save Outcome ----

# Save Plot ----

# Session Info ----
sessioninfo::session_info()
