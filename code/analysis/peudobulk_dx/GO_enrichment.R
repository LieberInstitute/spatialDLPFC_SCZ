# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(DOSE)
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

sig_gene <- gene_df |>
  filter(fdr_ntc <= 0.05) |>
  pull(ensembl)


query_string <-
paste0(sig_gene, collapse = "::")


dir.create(
  here("processed-data/PB_dx_genes/enrichment"),
  showWarnings = FALSE
)

# GO enrichment ----
## Overrepresentative analysis ----
ego <- enrichGO(
  gene = sig_gene,
  universe = gene_df$ensembl,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# dotplot(ego)

write.csv(
  ego@result,
  here(
    "processed-data/PB_dx_genes/enrichment",
    "GO_ERA_PRECAST_07.csv"
  )
)


## gene set enrichment analysis ----
geneList <- gene_df$logFC_scz
names(geneList) <- gene_df$ensembl
geneList <- sort(geneList, decreasing = TRUE)

go_gsea <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  minGSSize = 100,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  keyType = "ENSEMBL"
)

go_gsea@result$core_enrichment_symbol <- go_gsea@result$core_enrichment |> sapply(
  FUN = function(x) {
    # browser()
    str_split(x, "/")[[1]] |>
      bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") |>
      pull(SYMBOL) |>
      paste0(collapse = "/")
  }
)
#  setReadable(go_gsea, 'org.Hs.eg.db')


write.csv(
  go_gsea@result,
  here(
    "processed-data/PB_dx_genes/enrichment",
    "GO_GSEA_PRECAST_07.csv"
  )
)

# Disease Enrichment
## Over-representation analysis
library(DOSE)

sig_gene_entrez <- bitr(sig_gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

DO_ora <- enrichDO(
  gene = sig_gene_entrez$ENTREZID,
  universe = bitr(gene_df$ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID,
  ont = "DO",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.05,
  readable = TRUE,
)

DO_ora@result |>
  filter(p.adjust <= 0.05) |>
  write_csv(
    file = here(
      "processed-data/PB_dx_genes/enrichment",
      "DO_ERA_PRECAST_07.csv"
    )
  )

## gsea analysis
names(geneList) <- bitr(names(geneList),
 fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

DO_gsea <- gseDO(
  geneList = geneList,
  minGSSize = 120,
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  verbose = FALSE
)

setReadable(DO_gsea, 'org.Hs.eg.db')@result |>
write_csv(
    file = here(
      "processed-data/PB_dx_genes/enrichment",
      "DO_GSEA_PRECAST_07.csv"
    )
  )

# Session Info ----
session_info()
