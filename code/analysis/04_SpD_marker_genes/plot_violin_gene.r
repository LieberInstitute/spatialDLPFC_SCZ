# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(scater)
  library(sessioninfo)
})

# Load data -----
# load psueobulked data
## SpD07 ----
.spd <- "PRECAST_07"

sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)

# Select gene ----
## Oligodentrocytes ----
# from PEC data
oligo_genes <- c("PLP1", "MBP", "MOBP")
neuronal_genes <- c("SNAP25", "SLC17A7")

# stopifnot(length(oligo_genes) == length(oligo_ensembl))

# oligo_p <- plotExpression(
#   sce_pseudo,
#   features = Oligo_genes,
#   x = .spd,
#   color = "dx",
#   show_median = TRUE,
#   scales = "free_y",
#   swap_rownames = "gene_name"
# ) + scale_x_discrete(
#   limits = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4))
# ) +
#   scale_color_manual(
#     name = "Diagnosis",
#     values = c(
#       "ntc" = "blue", "scz" = "red"
#     ),
#     labels = c("NTC", "SCZ")
#   ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Neuronal genes ----


paneled_p <- plotExpression(
  sce_pseudo,
  features = c(neuronal_genes, Oligo_genes),
  x = .spd,
  color = "dx",
  show_median = TRUE,
  scales = "free_y",
  swap_rownames = "gene_name"
) + scale_x_discrete(
  limits = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4))
  # TODO: add the layer annotation here
) +
  scale_color_manual(
    name = "Diagnosis",
    values = c(
      "ntc" = "blue", "scz" = "red"
    ),
    labels = c("NTC", "SCZ")
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x.bottom = element_blank())


# Save plot -----


ggsave(
  filename = here("plots/04_SpD_marker_genes", "violin_plot_spd01-04.pdf"),
  plot = paneled_p,
  width = 12, height = 6
)




# Session Info ----
session_info()
