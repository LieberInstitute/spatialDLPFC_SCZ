# NOTE: deprecated code
#       For more info, check out `code/analysis/pseudobulk_dx/interaction_mdl/plot_spaghetti_.r`

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(limma)
  library(ggplot2)
  library(sessioninfo)
  library(SingleCellExperiment)
})

# Load Libray ----
## gene info ----
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
) |> select(
  ensembl, gene
)


sig_gene <- readxl::read_excel(
  here(
    "code/analysis/pseudobulk_dx",
    # "Test_90DEGs.xlsx"
    "Test_68DEGs.xlsx"
  ),
  col_names = FALSE
)[[1]]

# Log FC ----
cont.mat <- rbind(
  rep(-1, 7),
  rep(1, 7),
  matrix(0, nrow = 23, ncol = 7),
  cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
)
colnames(cont.mat) <- sprintf("spd%02d", 1:7)

fit <- readRDS(
  here(
    "processed-data/PB_dx_genes",
    "test_inter_PRECAST_07_20240627.rds"
  )
)

contrast_fit <- contrasts.fit(fit, cont.mat)
contrast_fit <- eBayes(contrast_fit)

cont_df <- topTable(contrast_fit, coef = sprintf("spd%02d", 1:7), num = Inf)

# Subset columns ----
fc_mat <- cont_df |>
  # dplyr::slice(c(1:10, (n() - 10):n())) |>
  select(starts_with("spd"))

centered_fc_mat <- fc_mat |>
  apply(1, scale, scale = FALSE) |>
  t()
colnames(centered_fc_mat) <- colnames(fc_mat)

# Assemble data matrix ----
raw_df <- cbind(
  fc_mat |> data.frame() |>
    rownames_to_column(var = "gene_id"),
  cont_df |> select(adj.P.Val)
) |>
  # Find gene name
  left_join(gene_df,
    by = c("gene_id" = "ensembl")
  )



turn_long_form <- function(dat) {
  dat |>
    data.frame() |>
    # rownames_to_column(var = "gene_id") |>
    pivot_longer(cols = starts_with("spd"))
}



# 69 genes ----
## Read Data ----
ful_df <- raw_df |> filter(gene %in% sig_gene)


bg_gene_df <- raw_df |>
  filter(gene %in% sig_gene) |>
  arrange(adj.P.Val) |>
  slice_tail(n = 10)


bg_gene_df |>
  turn_long_form() |>
  ggplot(aes(x = name, y = value, group = gene)) +
  # geom_point(aes(color = "non_sig"), alpha = 0.4) +
  geom_line(
    aes(color = "non_sig"),
    alpha = 0.4, linewidth = 0.5
  ) +
  geom_line(
    # data = ful_df |> filter(gene == "SST") |> turn_long_form(),
    data = ful_df |> filter(gene == "ALDH1A1") |> turn_long_form(),
    aes(color = "Sig (0.05)")
  ) +
  scale_x_discrete(
    limits = c(
      "spd07", "spd06", "spd02",
      "spd05", "spd03", "spd01", "spd04"
    )
  ) +
  labs(
    title = paste0("68 Genes"),
    y = "log2FC",
    x = ""
  ) +
  theme_light(base_size = 12)




# 172 genes





# Center across spd ----

# Divergent


# Session Info ----
session_info()
