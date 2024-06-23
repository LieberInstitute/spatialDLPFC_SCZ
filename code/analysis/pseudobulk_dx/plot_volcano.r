suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sessioninfo)
})

# Load Data ----

p_cut_off <- 0.10

gene_df <- read_csv(
  #   here(
  #     "processed-data/PB_dx_genes/",
  #     "test_PRECAST_07.csv"
  #   )
  # )

  "~/Downloads/test_PRECAST_07.csv"
) |> mutate(
  sig_05 = fdr_scz <= 0.05,
  sig_10 = fdr_scz <= 0.10
)




# sig_gene_df <- gene_df |>
#   filter(fdr_ntc <= p_cut_off)





# Make volcano plot ----
ggplot(
  gene_df,
  aes(x = logFC_scz, y = -log10(p_value_scz), color = sig_10)
) +
  geom_point() +
  # geom_label_repel(
  #   data = impl_gene_df, # Add labels last to appear as the top layer
  #   aes(label = gene),
  #   force = 2,
  #   nudge_y = 0.1
  # ) +
  # geom_label_repel(
  #   data = out_gene_df,
  #   aes(label = gene),
  #   color = "red",
  #   force = 2,
  #   nudge_y = -0.1
  # ) +
  labs(
    title = "Bulk analysis pooling all spots"
  ) +
  theme_minimal()


# Session Info ----
session_info()