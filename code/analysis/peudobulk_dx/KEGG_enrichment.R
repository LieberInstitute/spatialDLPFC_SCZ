# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  # library(DOSE)
  library(sessioninfo)
})

# Load Data ----
gene_df <- read_csv(
  #   here(
  #     "processed-data/PB_dx_genes/",
  #     "test_PRECAST_07.csv"
  #   )
  # )

  "~/Downloads/test_PRECAST_07.csv"
)




entrezid_df <- bitr(gene_df$ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

gene_df <- gene_df |> 
left_join(
  entrezid_df,
  by = c("ensembl" = "ENSEMBL")
)


sig_gene_df <- gene_df |>
  filter(fdr_ntc <= 0.05)

sig_gene <- sig_gene_df |>
  pull(ensembl)


# Export the genes for website
# query_string <- paste0(sig_gene, collapse = "::")


dir.create(
  here("processed-data/PB_dx_genes/enrichment"),
  showWarnings = FALSE
)


# KEGG & ClusterProfiler prep ----
search_kegg_organism("hsa", by = "kegg_code")

# Up-reg genes ----
up_gene <- gene_df |>
  filter(
    fdr_ntc <= 0.05,
    logFC_scz > 0
  ) |>
  pull(ENTREZID)

# TODO: fix this
KEGG_up_ora <- enrichKEGG(
  gene = up_gene,
  universe = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# TODO: make the results readable
# setReadable(KEGG_ora, 'org.Hs.eg.db')@result |>
KEGG_up_ora@result |>
write_csv(
    file = here(
      "processed-data/PB_dx_genes/enrichment",
      "KEGG_ORA_up_gene_PRECAST_07.csv"
    )
  )


# Down-reg genes ----
down_gene <- gene_df |>
  filter(
    fdr_ntc <= 0.05,
    logFC_scz < 0
  ) |>
  pull(ENTREZID)

length(down_gene)
# 61 up reg genes

KEGG_down_ora <- enrichKEGG(
  gene = down_gene,
  universe = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

KEGG_down_ora@result |>
write_csv(
    file = here(
      "processed-data/PB_dx_genes/enrichment",
      "KEGG_ORA_down_gene_PRECAST_07.csv"
    )
  )


# All genes ----
sig_gene_df <- gene_df |>
  filter(fdr_ntc <= 0.05)

sig_gene <- sig_gene_df |>
  pull(ENTREZID)


KEGG_all_ora <- enrichKEGG(
  gene = sig_gene,
  universe = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

KEGG_all_ora@result |>
write_csv(
    file = here(
      "processed-data/PB_dx_genes/enrichment",
      "KEGG_ORA_all_gene_PRECAST_07.csv"
    )
  )

# Session Info ----
session_info()
