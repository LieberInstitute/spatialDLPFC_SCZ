# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readxl)
  library(sessioninfo)
})

# Load Data ----
ruzicka_xlsx_path <- here(
  "code/analysis/12_cross-study_enrichment",
  "ruzicka_cell_type_DEG.xlsx"
)

## Load DEGs ----
ruzicka_sheets <- excel_sheets(ruzicka_xlsx_path)

# Error prevention
stopifnot(
  length(ruzicka_sheets) == 25
)

# format to list of data frames for each cell types
ruzicka_deg_list <- ruzicka_sheets |>
  set_names() |>
  map(
    ~ read_xlsx(
      ruzicka_xlsx_path,
      sheet = .x
    ) |>
      dplyr::select(gene, starts_with("Meta_")) |>
      mutate(
        cell_type = .x,
        ruzicka_sig_gene = case_when(
          Meta_adj.P.Val < 0.05 & abs(Meta_logFC) > 0.1 ~ TRUE,
          TRUE ~ FALSE
        ),
        # convert gene symbol to ensembl ID
        ensembl = AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = gene,
          column = "ENSEMBL",
          keytype = "SYMBOL",
          multiVals = "first"
        )
      )
  )

# Save data as RDS ----
saveRDS(
  ruzicka_deg_list,
  file = here(
    "processed-data/rds/12_cross-study_enrichment",
    "ruzicka_cell_type_DEG.rds"
  )
)

# Session info ----
session_info()
