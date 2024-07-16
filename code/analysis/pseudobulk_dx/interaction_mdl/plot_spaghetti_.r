# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(sessioninfo)
})

# Load data ----
## LogFC mat
log2fc_mat <- readRDS(
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_laminar_specific_log2FC.rds"
  )
)





# Make spaghetti plot -----

make_spaghetti_plot <- function(
    fc_mat,
    bg_genes,
    genes,
    .center = FALSE) {
  if (.center) {
    spd_comp <- fc_mat |>
      select(starts_with("spd"))

    .fc_mat <- spd_comp |>
      apply(1, scale, scale = FALSE) |>
      t()
    colnames(.fc_mat) <- colnames(spd_comp)

    fc_mat <- cbind(
      fc_mat |> select(-starts_with("spd")),
      .fc_mat
    )
  }

  turn_long_form <- function(dat) {
    dat |>
      data.frame() |>
      # rownames_to_column(var = "gene_id") |>
      pivot_longer(cols = starts_with("spd"))
  }

  fc_mat_long <- fc_mat |> turn_long_form()

  # Make plot
  ggplot(mapping = aes(x = name, y = value, group = gene)) +
    # Dashed line as reference
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.4) +
    # Create transparent paths for background genes

    geom_line(
      aes(color = "Background genes"),
      data = fc_mat_long |> filter(gene %in% bg_genes),
      alpha = 0.4, linewidth = 0.5
    ) +
    geom_point(
      data = fc_mat_long |> filter(gene %in% genes),
      aes(color = "Sig.")
    ) +
    geom_line(
      data = fc_mat_long |> filter(gene %in% genes),
      aes(color = "Sig.")
    ) +
    scale_x_discrete(
      limits = c(
        "spd07", "spd06", "spd02",
        "spd05", "spd03", "spd01", "spd04"
      )
    ) +
    scale_color_manual(
      values = c(
        "Background" = "gray",
        "Sig." = "red"
      )
    ) +
    labs(
      title = paste0("Spaghetti plot - ", genes),
      y = ifelse(.center, "log2FC (centered)", "log2FC"),
      x = ""
    ) +
    theme_light(base_size = 12)
}

treg_genes <- c("AKT3", "MALAT1", "ARID1B")

make_spaghetti_plot(
  fc_mat = log2fc_mat,
  bg_genes = treg_genes,
  genes = "ALDH1A1",
  .center = FALSE
)

make_spaghetti_plot(
  fc_mat = log2fc_mat,
  bg_genes = treg_genes,
  genes = "SST",
  .center = TRUE
)

## Define background genes





# gene info ----
# gene_df <- read_csv(
#   here(
#     "processed-data/PB_dx_genes/",
#     "test_PRECAST_07.csv"
#   )
# ) |> select(
#   ensembl, gene
# )


# sig_gene <- readxl::read_excel(
#   here(
#     "code/analysis/pseudobulk_dx",
#     # "Test_90DEGs.xlsx"
#     "Test_68DEGs.xlsx"
#   ),
#   col_names = FALSE
# )[[1]]

# # Make plot ----
# cont.mat <- rbind(
#   rep(-1, 7),
#   rep(1, 7),
#   matrix(0, nrow = 23, ncol = 7),
#   cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
# )
# colnames(cont.mat) <- sprintf("spd%02d", 1:7)

# fit <- readRDS(
#   here(
#     "processed-data/PB_dx_genes",
#     "test_inter_PRECAST_07_20240627.rds"
#   )
# )

# contrast_fit <- contrasts.fit(fit, cont.mat)
# contrast_fit <- eBayes(contrast_fit)

# cont_df <- topTable(contrast_fit, coef = sprintf("spd%02d", 1:7), num = Inf)

# # Subset columns ----
# fc_mat <- cont_df |>
#   # dplyr::slice(c(1:10, (n() - 10):n())) |>
#   select(starts_with("spd"))

# centered_fc_mat <- fc_mat |>
#   apply(1, scale, scale = FALSE) |>
#   t()
# colnames(centered_fc_mat) <- colnames(fc_mat)


# # Assemble data matrix ----
# raw_df <- cbind(
#   centered_fc_mat |> data.frame() |>
#     rownames_to_column(var = "gene_id"),
#   cont_df |> select(adj.P.Val)
# ) |>
#   # Find gene name
#   left_join(gene_df,
#     by = c("gene_id" = "ensembl")
#   )



# turn_long_form <- function(dat) {
#   dat |>
#     data.frame() |>
#     # rownames_to_column(var = "gene_id") |>
#     pivot_longer(cols = starts_with("spd"))
# }



# # 69 genes ----
# ## Read Data ----
# # ful_df <- raw_df |> filter(gene %in% c("ALDH1A1", "SOD2", bg_gene))


# bg_gene_df <- raw_df |> filter(gene %in% c(bg_gene))

# treg_gene_df <- raw_df |> filter(gene %in% c("AKT3", "MALAT1", "ARID1B"))


# treg_gene_df |>
#   turn_long_form() |>
#   ggplot(aes(x = name, y = value, group = gene)) +
#   # geom_point(aes(color = "non_sig"), alpha = 0.4) +
#   geom_line(
#     aes(color = "non_sig"),
#     alpha = 0.4, linewidth = 0.5
#   ) +
#   geom_line(
#     data = raw_df |> filter(gene %in% c("ALDH1A1", "SOD2")) |> turn_long_form(),
#     aes(color = "Sig (0.05)")
#   ) +
#   scale_x_discrete(
#     limits = c(
#       "spd07", "spd06", "spd02",
#       "spd05", "spd03", "spd01", "spd04"
#     )
#   ) +
#   labs(
#     title = "68 Genes",
#     y = "log2FC (Centered)",
#     x = ""
#   ) +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   theme_light(base_size = 12)


# ## Overlap with our genes
# sig_gene <- readxl::read_excel(
#   here(
#     "code/analysis/pseudobulk_dx",
#     # "Test_90DEGs.xlsx"
#     "Test_68DEGs.xlsx"
#   ),
#   col_names = FALSE
# )[[1]]


# intersect(res$gene, sig_gene)

# # Adjustment among 170 Sig_genes ----
# bg_gene <- res |>
#   filter(
#     gene %in% sig_gene
#   ) |>
#   mutate(sig_gene_adj_p = p.adjust(P.Value)) |>
#   arrange(sig_gene_adj_p, P.Value) |>
#   # View()
#   slice_tail(n = 10) |>
#   pull(gene)
