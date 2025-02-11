# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

spd_rds <- list.files(
  here(
    "processed-data", "rds", "layer_spd"
  ),
  pattern = ".rds"
)

# PCA analysis ----

dir.create(
  path = here(
    "plots/PB_DE", "SpD_PB_PCA"
  ),
  showWarnings = FALSE
)

for (.file in spd_rds) {
  # .file <- spd_rds[1]

  .spd <- str_remove(.file, "test_spe_pseudo_") |>
    str_remove(".rds")

  sce_pseudo <- readRDS(
    here(
      "processed-data", "rds", "layer_spd",
      .file
    )
  )


  set.seed(20240411)
  sce_pseudo <- runPCA(sce_pseudo)

  pdf(
    here(
      "plots/PB_DE", "SpD_PB_PCA",
      paste0("test_", .spd, ".pdf")
    )
  )
  sce_pseudo <- runPCA(sce_pseudo)
  for (.var in c("age", "sex", "spd", "dx", "slide_id", "lot_num")) {
    plotPCA(
      sce_pseudo,
      colour_by = .var,
      ncomponents = 6,
      point_size = 0.3,
      label_format = c("%s %02i", " (%i%%)")
    ) |> print()
  }

  plotExplanatoryVariables(sce_pseudo,
    variables = c("dx", "age", "sex", "slide_id", "spd", "sample_id")
  )
  dev.off()
}



# # Dx DE analysis ----
# pdf(
#   here(
#     "plots/PB_DE",
#     "test_dx_volcano_plot.pdf"
#   )
# )

# for (.file in spd_rds) {
#   # .file <- spd_rds[1]

#   .spd <- str_remove(.file, "test_spe_pseudo_") |>
#     str_remove(".rds")

#   sce_pseudo <- readRDS(
#     here(
#       "processed-data", "rds", "layer_spd",
#       .file
#     )
#   )
#   dx_mod <-
#     registration_model(
#       sce_pseudo,
#       covars = c("spd", "age", "sex", "slide_id"),
#       var_registration = "dx"
#     )

#   dx_block_cor <-
#     registration_block_cor(
#       sce_pseudo,
#       registration_model = dx_mod,
#       var_sample_id = "sample_id"
#     )

#   dx_res <- registration_stats_enrichment(
#     sce_pseudo,
#     block_cor = dx_block_cor,
#     covars = c("spd", "age", "sex", "slide_id"),
#     var_registration = "dx",
#     gene_ensembl = "gene_id",
#     gene_name = "gene_name"
#   )

#   dir.create(
#     here("processed-data/PB_dx_genes"),
#     showWarnings = FALSE
#   )

#   ## Save DE results ----
#   dx_res |>
#     arrange(fdr_scz) |>
#     write.csv(
#       here(
#         "processed-data/PB_dx_genes",
#         paste0("test_", .spd, ".csv")
#       )
#     )

#   ## Volcano Plot ----
#   out_gene_df <- dx_res |>
#     arrange(fdr_scz) |>
#     slice_head(n = 1)

#   impl_gene_df <- dx_res |>
#     filter(gene %in% c(
#       "PVALB",
#       "NOS1",
#       "SST",
#       "CHODL",
#       "GRIN2A",
#       "SV2A",
#       "DLG4",
#       "C4A",
#       "C3"
#     )) |>
#     select(ensembl, gene, ends_with("scz"))

#   n_sig_gene <- dx_res |>
#     filter(fdr_scz <= 0.05) |>
#     nrow()

#   print(
#     ggplot(
#       dx_res,
#       aes(
#         x = logFC_scz, y = -log10(fdr_scz),
#         color = fdr_scz <= 0.05
#       )
#     ) +
#       geom_point(alpha = 0.8) +
#       geom_label_repel(
#         data = impl_gene_df,
#         aes(label = gene),
#         force = 2,
#         nudge_y = 0.1
#       ) +
#       geom_label_repel(
#         data = out_gene_df,
#         aes(label = gene),
#         force = 2,
#         nudge_y = -0.1
#       ) +
#       labs(
#         title = paste0(
#           "Pseudobulk analysis by Dx - ", .spd,
#           " ( ", n_sig_gene, " sig genes)"
#         )
#       ) +
#       theme_minimal()
#   )
# }
# dev.off()

# Session Info ----
sessioninfo::session_info()
