# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggplot2)
  library(sessioninfo)
})

# Load data ----
## Synthetic null
null_mat <- readRDS(
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_permed_laminar_specific_log2FC.rds"
  )
) |> slice_tail(n = 50)

## LogFC mat
log2fc_mat <- readRDS(
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_laminar_specific_log2FC.rds"
  )
)

gene_172_ensembl <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
) |>
  filter(fdr_scz <= 0.10) |>
  pull(ensembl)

log2fc_mat_172 <- log2fc_mat |> filter(ensembl %in% gene_172_ensembl)




# Make spaghetti plot -----
make_spaghetti_plot_null <- function(
    fc_mat,
    bg_mat,
    bg_genes,
    genes,
    .center = FALSE) {
  fc_mat <- rbind(
    fc_mat,
    bg_mat
  )

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

  ret_p <- ggplot(mapping = aes(x = name, y = value, group = gene)) +
    # Dashed line as reference
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.4) +
    # Create transparent paths for background genes
    geom_line(
      aes(color = "Synthesized"),
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
        "Synthesized" = "gray",
        "Sig." = "red"
      )
    ) +
    labs(
      title = paste0("Spaghetti plot - ", genes),
      y = ifelse(.center, "log2FC (centered)", "log2FC"),
      x = ""
    ) +
    theme_light(base_size = 12)

  return(ret_p)
}

# treg_genes <- c("AKT3", "MALAT1", "ARID1B")




pdf(
  here(
    "plots/PB_dx_genes/spaghetti_plot",
    "spatially_divergent_genes_syn_background.pdf"
  ),
  height = 5
)
make_spaghetti_plot_null(
  fc_mat = log2fc_mat_172,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "HNRNPH3",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat_172,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "ALDH1A1",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat_172,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "SOD2",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat_172,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "CCNI",
  .center = FALSE
) |> print()
dev.off()


# Debug why some permuted genes still have very large logFC
# ENSG00000228543


make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "AC003684.1",
  .center = FALSE
)

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "TOM1L1",
  .center = FALSE
)


# Housekeeping genes ----
# Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2950262/
# Note: "PBDG" are not in the data set


pdf(
  here(
    "plots/PB_dx_genes/spaghetti_plot",
    "spatially_divergent_genes_showing House Keeping gene.pdf"
  ),
  height = 5
)
make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "GUSB",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "PPIA",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "UBC",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "SDHA",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "ACTB",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "TBP",
  .center = FALSE
) |> print()

make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "B2M",
  .center = FALSE
) |> print()


make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = null_mat$gene,
  genes = "GAPDH",
  .center = FALSE
) |> print()
dev.off()





pdf(
  here(
    "plots/PB_dx_genes/spaghetti_plot",
    "spatially_divergent_genes_HK_background.pdf"
  ),
  height = 5
)
print(
  make_spaghetti_plot_null(
    fc_mat = log2fc_mat,
    bg_mat = null_mat,
    bg_genes = c("GUSB", "PPIA", "UBC", "SDHA", "ACTB", "TBP", "B2M", "GAPDH"),
    genes = "HNRNPH3",
    .center = FALSE
  ) +
    scale_color_manual(
      values = c(
        "Synthesized" = "gray",
        "Sig." = "red"
      ),
      labels = c(
        "Synthesized" = "Housekeeper",
        "Sig." = "Sig."
      )
    )
)

print(
  make_spaghetti_plot_null(
    fc_mat = log2fc_mat,
    bg_mat = null_mat,
    bg_genes = c("GUSB", "PPIA", "UBC", "SDHA", "ACTB", "TBP", "B2M", "GAPDH"),
    genes = "ALDH1A1",
    .center = FALSE
  ) +
    scale_color_manual(
      values = c(
        "Synthesized" = "gray",
        "Sig." = "red"
      ),
      labels = c(
        "Synthesized" = "Housekeeper",
        "Sig." = "Sig."
      )
    )
)

print(
  make_spaghetti_plot_null(
    fc_mat = log2fc_mat,
    bg_mat = null_mat,
    bg_genes = c("GUSB", "PPIA", "UBC", "SDHA", "ACTB", "TBP", "B2M", "GAPDH"),
    genes = "SOD2",
    .center = FALSE
  ) +
    scale_color_manual(
      values = c(
        "Synthesized" = "gray",
        "Sig." = "red"
      ),
      labels = c(
        "Synthesized" = "Housekeeper",
        "Sig." = "Sig."
      )
    )
)

print(make_spaghetti_plot_null(
  fc_mat = log2fc_mat,
  bg_mat = null_mat,
  bg_genes = c("GUSB", "PPIA", "UBC", "SDHA", "ACTB", "TBP", "B2M", "GAPDH"),
  genes = "CCNI",
  .center = FALSE
) +
  scale_color_manual(
    values = c(
      "Synthesized" = "gray",
      "Sig." = "red"
    ),
    labels = c(
      "Synthesized" = "Housekeeper",
      "Sig." = "Sig."
    )
  ))
dev.off()


# Session Info ----
sessioninfo::session_info()
