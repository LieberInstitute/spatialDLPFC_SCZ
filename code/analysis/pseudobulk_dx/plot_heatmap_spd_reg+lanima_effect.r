library(here)


# Load Spd Registration
source(
  here(
    "code/analysis/pseudobulk_dx",
    "plot_sig_gene_spd_enrichment.r"
  )
)


# Load Effect Size plot
source(
  here(
    "code/analysis/pseudobulk_dx/interaction_mdl",
    "plot_heatmap_effect_size.R"
  )
)


# Check if the order of the two heatmap is the same.
stopifnot(
  all(
    heatmap_all@row_names_param$labels ==
      effect_heatmap@row_names_param$labels
  )
)

# stopifnot(
#   all(
#     heatmap_all@row_names_param$labels[heatmap_all |>
#       row_order() |>
#       unlist() |>
#       c()] ==
#       effect_heatmap@row_names_param$labels[effect_heatmap |>
#         row_order() |>
#         unlist()]
#   )
# )


pdf(
  here(
    "plots/PB_dx_genes",
    "test_Heatmap_spd_reg_N_effect_size.pdf"
  ),
  height = 20
)
heatmap_all + effect_heatmap
dev.off()
