suppressPackageStartupMessages({
  library(tidyverse)
  library(STRINGdb)
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

## Format result to STRINGdb example data ----
# data(diff_exp_example1)
# diff_exp_example1 |> str()
# r$> diff_exp_example1 |> str()
# 'data.frame':   20861 obs. of  3 variables:
#  $ pvalue: num  0.000102 0.000139 0.000172 0.000174 0.000199 ...
#  $ logFC : num  3.33 3.82 3.31 3.02 3.85 ...
#  $ gene  : chr  "VSTM2L" "TBC1D2" "LENG9" "TMEM27" ...

diff_exp_df <- gene_df |>
  select(
    pvalue = fdr_scz, logFC = logFC_scz,
    gene, ensembl
  ) |>
  arrange(pvalue) |>
  data.frame() # STRINGdb doesn't work with tibble




# STRING db Enrichment ---
## Laod STRING data base ----
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606, # Human
  score_threshold = 400,
  network_type = "full",
  input_directory = ""
)


## Analysis using default setting ----
# deg_mapped_gene <- string_db$map(
#   diff_exp_df,
#   my_data_frame_id_col_names = "gene",
#   removeUnmappedRows = TRUE
# )
# # Warning:  we couldn't map to STRING 8% of your identifiers

# # deg_mapped_ensembl <-  string_db$map(
# #   diff_exp_df,
# #   my_data_frame_id_col_names = "ensembl",
# #   removeUnmappedRows = TRUE
# # )
# # There's 8% gene can't be mapped when using ensembleID

# hits <- deg_mapped_gene$STRING_id[1:200]
# string_db$plot_network(hits)


## Analysis with significant subset of genes ---
sig_gene_df <- diff_exp_df |>
  filter(pvalue <= 0.05) # STRINGdb doesn't work with tibble

sig_gene_df |>
  string_db$map(
    my_data_frame_id_col_names = "gene",
    removeUnmappedRows = TRUE
  ) |>
  pull(STRING_id) |>
  string_db$get_link()

sig_mapped_gene <- string_db$map(
  sig_exp_df,
  my_data_frame_id_col_names = "gene",
  removeUnmappedRows = TRUE
)


hits <- sig_mapped_gene$STRING_id
string_db$plot_network(hits)

### Clustering (Deprecated) ----
# clustersList <- string_db$get_clusters(hits)
# length(clustersList)
# par(mfrow = c(4, 4))
# for (i in seq_along(clustersList)) {
#   string_db$plot_network(clustersList[[2]])
# }
# par(mfrow = c(1, 1))

# head(diff_exp_example1)


## Up Gene ----
up_gene_df <- sig_gene_df |> filter(logFC > 0)
up_mapped_gene <- string_db$map(
  up_gene_df,
  my_data_frame_id_col_names = "gene",
  removeUnmappedRows = TRUE
)

string_db$get_link(up_mapped_gene$STRING_id)
# string_db$plot_network(up_mapped_gene$STRING_id)


## Down Gene ----
down_mapped_gene <- sig_gene_df |>
  filter(logFC < 0) |>
  string_db$map(
    my_data_frame_id_col_names = "gene",
    removeUnmappedRows = TRUE
  )
# string_db$plot_network(down_mapped_gene$STRING_id)

string_db$get_link(down_mapped_gene$STRING_id)

# Session Info ----
session_info()
