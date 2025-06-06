# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(AnnotationDbi)
  # library(org.Hs.eg.db)
  library(readxl)
  library(sessioninfo)
  library(biomaRt)
})

# Load Data ----


## Load Ensembl Mart ----
# Connect to the Ensembl database
ensembl <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
) # for human genes

## Load DEGs ----
ruzicka_xlsx_path <- here(
  "code/analysis/12_cross-study_enrichment",
  "ruzicka_cell_type_DEG.xlsx"
)

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
        )
      )
  )


# Idenitfy ensembl id for all gene symbols ----
all_gene_symbols <- ruzicka_deg_list |>
  map(
    ~ .x |> pull(gene)
  ) |>
  unlist() |>
  unique()

result <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = all_gene_symbols,
  mart = ensembl
)

# Check if all gene symbols were mapped to Ensembl IDs
stopifnot(sum(is.na(result$ensembl_gene_id)) == 0)


# Merge Ensembl IDs with ruzicka_deg_list ----
fnl_ruzicka_deg_list <- ruzicka_deg_list |>
  map(
    ~ .x |>
      left_join(
        result,
        by = c("gene" = "hgnc_symbol")
      ) |>
      dplyr::rename(
        ensembl = ensembl_gene_id
      )
  )


# Save data as RDS ----
saveRDS(
  fnl_ruzicka_deg_list,
  file = here(
    "processed-data/rds/12_cross-study_enrichment",
    "ruzicka_cell_type_DEG.rds"
  )
)

# Session info ----
session_info()
