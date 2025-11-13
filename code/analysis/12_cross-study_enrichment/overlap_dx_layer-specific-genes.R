# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load Data ----
## 172 genes ----
spd_deg_172 <- read_csv(
  here(
    "code/analysis/10_dx_deg_adjust_spd",
    "172_prelim_fdr010.csv"
  ),
  col_types = cols()
)

## Load Layer-specific DEGs ----
spd_files <- list.files(
  here("processed-data/rds/11_dx_deg_interaction"),
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

names(spd_files) <- str_extract(
  spd_files,
  "(?<=layer_specific_logFC_).*?(?=\\.csv)"
)

spd_deg_df <- imap_dfr(
  spd_files,
  ~ read_csv(.x, col_types = cols()) |>
    mutate(
      gene = AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = gene_id,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      ),
      spd = .y
    ) |>
    filter(
      P.Value < 0.05
    )
)

spd_deg_df <- spd_deg_df |>
  mutate(
    overlap_172 = ifelse(
      gene_id %in% spd_deg_172$ensembl,
      "in 172 genes",
      "Not in 172 genes"
    )
  )

write_csv(
  spd_deg_df,
  here(
    "code/analysis/12_deg_integration",
    "overlap_layer_specific_overall_genes.csv"
  )
)

spd_deg_df |> filter(overlap_172 == "in 172 genes") |> pull(gene) |> unique() |> length()


## with filter down genes ----

tmp <- read_csv(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_specific_genes_fdr.csv"
  )
) |>
  mutate(
    overlap_172 = ifelse(
      ENSEMBL %in% spd_deg_172$ensembl,
      "in 172 genes",
      "Not in 172 genes"
    )
  )

tmp |> filter(overlap_172 == "in 172 genes") |> pull(ENSEMBL) |> unique() |> length()

write_csv(
  tmp,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "layer_specific_genes_fdr_with_172_annotation.csv"
  )
)


# Session Info ----
sessioninfo::session_info()
