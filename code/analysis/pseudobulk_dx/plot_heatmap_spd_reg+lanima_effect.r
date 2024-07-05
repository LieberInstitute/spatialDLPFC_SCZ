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
  identical(
    heatmap_all@row_names_param$labels,
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

# CREATE dot plot for GO Pathway
source(
  here(
    "code/analysis/dx_deg_enrichment",
    "plot_enrich_diff_string.r"
  )
)


# Check if the order of the two heatmap is the same.
stopifnot(
  all(
    identical(
      heatmap_all@row_names_param$labels,
      heatmap_go@row_names_param$labels
    )
  )
)

# pdf(
#   here(
#     "plots/PB_dx_genes",
#     sprintf(
#       "test_Heatmap_spd_reg_N_effect_size_%02dGene.pdf",
#       n_gene
#     )
#   ),
#   height = 20
# )
# print(heatmap_all + effect_heatmap)
# dev.off()




pdf(
  here(
    "plots/PB_dx_genes",
    sprintf(
      "test_Heatmap_spd_reg_N_effect_size_%02dGene.pdf",
      n_gene
    )
  ),
  height = 20
)
print(heatmap_all + effect_heatmap + heatmap_go)

dev.off()
