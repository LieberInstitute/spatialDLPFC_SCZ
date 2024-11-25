# Load library ----
library(data.table)
library(fgsea)
library(ggplot2)

# Load Data ----
data(examplePathways)
data(exampleRanks)


str(examplePathways)
str(exampleRanks)

# fgsea analysis ----
fgseaRes <- fgsea(
  pathways = examplePathways,
  stats = exampleRanks,
  minSize = 15,
  maxSize = 500
)

# Viz ----
## Viz one pathway ----
plotEnrichment(
  examplePathways[["5991130_Programmed_Cell_Death"]],
  exampleRanks
) + labs(title = "Programmed Cell Death")

## Viz up and down pathway ----
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes,
  gseaParam = 0.5
)

# Session Info ----
sessioninfo::session_info()
