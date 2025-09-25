# Load library -----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(readxl)
  library(CellChat)
  library(circlize)
  library(ComplexHeatmap)
  library(sessioninfo)
})

# Load data -----
## Load SCZ-risk genes ----
# NOTE: supp table 12 from Trubetskoy
trubetskoy_df <- read_excel(
  here(
    "code/analysis/Layer_layer_communication/genetic_risk",
    "Supplementary Table 12.xlsx"
  ),
  sheet = "Prioritised"
)

symbols_gene <- trubetskoy_df |>
  # filter(gene_biotype == "protein_coding") |>
  pull(Symbol.ID)

## Load dx-DEG genes -----
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
)

sig_genes <- gene_df |> filter(fdr_scz <= 0.10)
stopifnot(length(sig_genes) != 172)


# Finding overlapping with dx DEG ----
intersect(sig_genes$ensembl, trubetskoy_df$Ensembl.ID)
intersect(sig_genes$gene, symbols_gene)

## Load CellChat Pathway -----
cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "merged_cellchat_annotated.rds"
  )
)

# CellChat Characterization -----
## Identify possible pathways in DB -----
### risk genes serves as both LR roles ----
LR_trubetskoy_both <- identifyOverExpressedInteractions(
  cellchat, # This can be any cellchat object that have @DB
  features = symbols_gene,
  return.object = FALSE, variable.both = TRUE
)

write_csv(
  LR_trubetskoy_both,
  here(
    "code/analysis/Layer_layer_communication/genetic_risk",
    "LR_trubetskoy_both.csv"
  )
)

### risk genes serves as either LR roles ----
LR_trubetskoy_single <- identifyOverExpressedInteractions(
  cellchat,
  features = symbols_gene,
  return.object = FALSE, variable.both = FALSE
)

write_csv(
  LR_trubetskoy_single,
  here(
    "code/analysis/Layer_layer_communication/genetic_risk",
    "LR_trubetskoy_single.csv"
  )
)

# quick statistics
nrow(LR_trubetskoy_single)
unique(LR_trubetskoy_single$pathway_name) |> paste(collapse = ",")
unique(LR_trubetskoy_single$pathway_name) |> length()

## Idnetify CellChat identified pathways ----
cellchat_pathways <- rankNet(cellchat, return.data = TRUE)$signaling.contribution



cellchat_diff_pathways <- unique(cellchat_pathways$name) |> as.character()

cellchat_diff_pathways |> length()

overlap_pathways <- intersect(
  cellchat_diff_pathways,
  LR_trubetskoy_single$pathway_name
)

overlap_pathways |> paste(collapse = ", ")

# Plot Pathways in NTC/SCZ and diff ----
## Load cellchat object of each cohort ----
ntc_cellchat <- readRDS(here(
  "processed-data/layer_layer_comm",
  "ntc_cellchat_annotated.rds"
))
scz_cellchat <- readRDS(
  here(
    "processed-data/layer_layer_comm",
    "scz_cellchat_annotated.rds"
  )
)

object.list <- list(ntc = ntc_cellchat, scz = scz_cellchat)

# Source code to create diff circular plot
source(here(
  "code/analysis/SPG_analysis/ccc_neun_claudin",
  "3_A-Sender_Receiver_Chord.r"
))
## Heatmap show diff
source(
  here(
    "code/analysis/Layer_layer_communication/vulnerablility_charact",
    "fun_extract_pathway_diff_mat.r"
  )
)

plot_pathway_panel <- function(vec_pathway, plot_path) {
  # TODO: create pdf file to save
  pdf(
    plot_path,
    height = 7, width = 9
  )
  for (.path in vec_pathway) {
    ## 2 circular plots for NTC and SCZ respectively

    par(mfrow = c(1, 2), xpd = TRUE)
    pathways.show <- .path
    weight.max <- getMaxWeight(object.list,
      slot.name = c("netP"),
      attribute = pathways.show
    ) # control the edge weig
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]],
        signaling = pathways.show, layout = "chord",
        edge.weight.max = weight.max[1], edge.width.max = 10,
        signaling.name = paste(pathways.show, names(object.list)[i])
      )
    }
    par(mfrow = c(1, 1), xpd = TRUE)
    ## Circular plots for diff
    my_netVisual_aggregate(ntc_cellchat, scz_cellchat,
      signaling = .path,
      measure = "weight", layout = "chord"
    )

    # Extract the pathway matrix
    pathway_matrix <- extract_diff_pathway(
      ntc_cellchat, scz_cellchat,
      signaling = .path, measure = "weight"
    ) |>
      as.matrix()

    # Create row and column annotations
    row_annotation <- rowAnnotation(
      hist = anno_barplot(pathway_matrix |> rowSums(), which = "row")
    )
    col_annotation <- HeatmapAnnotation(
      hist = anno_barplot(pathway_matrix |> colSums(), which = "column")
    )


    # Visualize the pathway matrix using a heatmap with annotations
    Heatmap(
      pathway_matrix,
      name = .path,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title = "Sender",
      column_title = "Receiver",
      col = colorRamp2(c(min(pathway_matrix), 0, max(pathway_matrix)), c("blue", "white", "red")),
      top_annotation = col_annotation,
      left_annotation = row_annotation
    ) |> print()
  }

  dev.off()
}

plot_pathway_panel(
  vec_pathway = overlap_pathways,
  plot_path = here(
    "code/analysis/Layer_layer_communication/genetic_risk",
    "trubetskoy_overlap_pathways_w_diff.pdf"
  )
)




# Session Info ----
session_info()
