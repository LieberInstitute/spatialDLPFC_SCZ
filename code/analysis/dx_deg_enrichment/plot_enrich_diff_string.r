library(here)
library(readr)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

intrst_terms <- c(
  "Response to stimulus",
  "Response to stress",
  "Cellular respiration",
  "Detoxification",
  "Immune system process",
  "Somatodendritic compartment",
  "Mixed, incl. Complement and coagulation cascades, and Protein-lipid complex"
)

sig_gene <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    # "Test_90DEGs.xlsx"
    "Test_68DEGs.xlsx"
  ),
  col_names = FALSE
)[[1]]

string_df <- readr::read_tsv(
  here("code/analysis/dx_deg_enrichment/string_functional_annotations.tsv")
) |> filter(
  `term description` %in% intrst_terms,
  `#node` %in% sig_gene
)

tmp_df <- string_df |>
  group_by(`#node`) |>
  summarize(terms = list(`term description`))

path_mat <- tmp_df$terms |>
  set_names(tmp_df$`#node`) |>
  list_to_matrix() |>
  t()


miss_genes <- setdiff(sig_gene, rownames(path_mat))

# neg_gene <- c("MALAT1", "ARID1B", "AKT3")
neg_gene <- NULL

miss_genes <- c(miss_genes, neg_gene)


miss_genes_mat <- matrix(
  0,
  nrow = length(miss_genes),
  ncol = length(intrst_terms)
)
rownames(miss_genes_mat) <- miss_genes

path_mat <- path_mat |> rbind(miss_genes_mat)



# gene_names_hc_ordered <- readRDS(
#   here(
#     "code/analysis/pseudobulk_dx",
#     sprintf(
#       "spd_hierarchical_cluster_order_%02d_gene.rds",
#       n_gene
#     )
#   )
# )

path_ht_mat <- path_mat[
  # gene_names_hc_ordered,
  c(sig_gene, neg_gene),
  intrst_terms
]

colnames(path_ht_mat)[length(intrst_terms)] <- "Mixed, complement & coagulation"


heatmap_go <- Heatmap(
  path_ht_mat,
  rect_gp = gpar(type = "none"),
  column_title = "Func. Annot.",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  # width = unit(30, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x = x, y = y, width = width, height = height,
      gp = gpar(col = "grey", fill = NA)
    )
    if (path_ht_mat[i, j]) {
      grid.circle(
        x = x, y = y,
        r = unit(1, "mm")
        #  r = min(unit.c(width, height))
        ,
        gp = gpar(fill = "black", col = NA)
      )
    }
  },
  show_column_names = TRUE, # Turn off legend
  show_row_names = TRUE, # Turn off legend,
  show_heatmap_legend = FALSE
)






#   col = colorRampPalette(c("black", "white"))(2),

# )
