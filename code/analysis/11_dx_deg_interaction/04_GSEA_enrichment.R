# Load required packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  # library(DOSE)
})


# Load Data ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

spd_deg_list <- map(
  spd_files,
  ~ read_csv(.x, col_types = cols())
)

names(spd_deg_list) <- str_extract(
  spd_files,
  "(?<=layer_specific_logFC_).*?(?=\\.csv)"
)



# gsea analysis ----
tmp_res <- spd_deg_list |>
  imap(.f = function(.df, idx) {
    # browser()
    # format gene list
    geneList <- .df |>
      arrange(desc(t)) |>
      dplyr::select(gene_id, t) |>
      deframe()

    set.seed(20250422)
    # Run GSEA
    ego_gsea <- gseGO(
      geneList = geneList,
      OrgDb = org.Hs.eg.db,
      ont = "BP", # focus on BP
      # ont = "ALL", # to have all
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      keyType = "ENSEMBL",
      pAdjustMethod = "fdr"
    )

    # make gene readable
    ego_gsea@result <- ego_gsea@result |>
      mutate(
        core_enrichment_symbol = core_enrichment |> sapply(
          FUN = function(x) {
            # browser()
            str_split(x, "/")[[1]] |>
              bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") |>
              pull(SYMBOL) |>
              paste0(collapse = "/")
          }
        )
      )

    # save results
    write_csv(
      ego_gsea@result,
      here(
        "processed-data/rds/11_dx_deg_interaction",
        paste0("gsea_int_", gsub("[^[:alnum:]_]", "_", idx), ".csv")
      )
    )
  })
