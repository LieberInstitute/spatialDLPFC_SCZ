# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(spatialLIBD)
  library(readxl)
  library(sessioninfo)
})

# Load data ----
## Load Ruzicka cell type DEGs ----
# List all sheets in the Ruzicka Excel file
ruzicka_sheets <- excel_sheets(
  here(
    "code/analysis/12_deg_integration",
    "ruzicka_cell_type_DEG.xlsx"
  )
)

# Error prevention
stopifnot(
  length(ruzicka_sheets) == 25
)

# format to list of data frames for each cell types
ruzicka_deg_list <- ruzicka_sheets |>
  set_names() |>
  map(
    ~ read_xlsx(
      here(
        "code/analysis/12_deg_integration",
        "ruzicka_cell_type_DEG.xlsx"
      ),
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

# Num of sig gene by cell type
ruzicka_deg_list |>
  map_int(
    ~ sum(.x$ruzicka_sig_gene, na.rm = TRUE)
  )


ruzicka_deg_df <- ruzicka_deg_list |>
  imap_dfr(~ .x |>
    filter(ruzicka_sig_gene) |>
    select(gene) |>
    mutate(cell_type = .y)) |>
  select(term = cell_type, gene)


## DLPFC data -----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  # pattern = "layer_specific_logFC_.*\\.csv",
  pattern = "layer_restricted_logFC_.*\\.csv",
  full.names = TRUE
)

names(spd_files) <- str_extract(
  spd_files,
  # "(?<=layer_specific_logFC_).*?(?=\\.csv)"
  "(?<=layer_restricted_logFC_).*?(?=\\.csv)"
)

spd_deg_list <- map(
  spd_files,
  ~ read_csv(.x, col_types = cols()) |>
    mutate(gene = AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ))
)


# enrich_res <- spd_deg_list |> map(
#   .f = function(.x) {
#     tmp <- .x |>
#       filter(P.Value < 0.05) |>
#       pull(gene)
#     names(tmp) <- NULL
#     clusterProfiler::enricher(
#       gene = tmp,
#       universe = spd_deg_list[[1]] |> pull(gene),
#       TERM2GENE = ruzicka_deg_df,
#       pAdjustMethod = "fdr",
#       pvalueCutoff = 1,
#       minGSSize = 0
#     )
#   }
# )

tmp <- spd_deg_list[[1]] |> filter(adj.P.Val < 0.10) |> pull(gene)
names(tmp) <- NULL


# Enrichment Test ----
clusterProfiler::enricher(
  gene = tmp,
  TERM2GENE = ruzicka_deg_df,
  # universe = spd_deg_list[[1]] |> pull(gene_id),
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.20
)

# Enrichment Test ----
# NOTE: using the universe to specific background genes doesn't change the number of bgRatio
clusterProfiler::enricher(
  gene = tmp,
  TERM2GENE = ruzicka_deg_df,
  universe = spd_deg_list[[1]] |> pull(gene_id),
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.20
)



# # Enrichment Test ----
clusterProfiler::enricher(
  gene = tmp,
  TERM2GENE = ruzicka_deg_df |>
    filter(!term %in% c("Ex-L2")),
  # universe = spd_deg_list[[1]] |> pull(gene_id),
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.20
)

ruzicka_deg_df |> group_by(term) |>
  summarise(
    n_genes = n(),
    n_sig_genes = sum(gene %in% tmp)
  )

unique(ruzicka_deg_df$gene)

ruzicka_deg_df |>
    filter(!term %in% c("Ex-L2")) |>
    pull(gene) |>
    unique() |>
    length()




ruzicka_deg_df$gene |>
  unique() |>
  length()




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


# Session Info ----
sessioninfo::session_info()
