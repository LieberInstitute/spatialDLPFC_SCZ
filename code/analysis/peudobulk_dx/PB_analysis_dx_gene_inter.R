# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(limma)
  # library(spatialLIBD)
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

# Dx DE analysis ----
pdf(
  here(
    "plots/PB_DE",
    "test_dx_volcano_plot.pdf"
  )
)
.file <- "test_spe_pseudo_PRECAST_07.rds"
{
# for (.file in spd_rds) {
  # .file <- spd_rds[1]

  .spd <- str_remove(.file, "test_spe_pseudo_") |>
    str_remove(".rds")

  sce_pseudo <- readRDS(
    here(
      "processed-data", "rds", "layer_spd",
      .file
    )
  )

  sce_pseudo$spd <- sce_pseudo[[.spd]]


  dx_mod <- model.matrix(~ dx*spd + age + sex + slide_id,
   colData(sce_pseudo))

  int_terms <- grep(":", dx_mod |> colnames(), value = TRUE)



  # Naive without adjusting for sample-id random effect
  fit <- lmFit(logcounts(sce_pseudo),
  design = dx_mod)
  fit <- eBayes(fit)

# Omnibus test for 
topTable(fit, coef = c("dxscz", int_terms), num = 30)

# Omnibus test for interaction terms
topTable(fit, coef = c(int_terms), num = 30)

# test of the deviation (FALSE)
topTable(fit, coef = int_terms[6]#, n=Inf
)

# Test for diagnosis
topTable(fit, coef = "dxscz", n=30
)


  # TODO: add one that adjusts effect


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
# }
# dev.off()

# Session Info ----
sessioninfo::session_info()
