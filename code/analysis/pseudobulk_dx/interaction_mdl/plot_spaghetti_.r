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

  ret_p <- ggplot(mapping = aes(x = name, y = value, group = gene)) +
    # Dashed line as reference
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.4) +
    # Create transparent paths for background genes
    # geom_line(
    #   aes(color = "Background genes"),
    #   data = fc_mat_long |> filter(gene %in% bg_genes),
    #   alpha = 0.4, linewidth = 0.5
    # ) +
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

  return(ret_p)
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

## Define background genes ----
# Laod the 172 sig Dx DEG at FDR 0.1
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
) |> mutate(
  sig_05 = fdr_scz <= 0.05,
  sig_10 = fdr_scz <= 0.10
)

dx_deg <- gene_df |>
  filter(fdr_scz <= 0.10) |>
  arrange(fdr_scz) |>
  pull(gene)

# Choose the genes ----

pdf(
  here("plots/PB_dx_genes/spaghetti_plot/",
  "spatially_divergent_genes_w_172bg.pdf")
)
ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "SST",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "SST",
    .center = FALSE
  ),
  nrow = 2
) |> print()

ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "FKBP5",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "FKBP5",
    .center = FALSE
  ),
  nrow = 2
) |> print()

ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "SOD2",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "SOD2",
    .center = FALSE
  ),
  nrow = 2
) |> print()

ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "ALDH1A1",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = tail(dx_deg, 20),
    genes = "ALDH1A1",
    .center = FALSE
  ),
  nrow = 2
) |> print()
dev.off()



bg_gene_all <- gene_df |> arrange(fdr_scz) |> 
slice_tail(n = 20) |> pull(gene)


pdf(
  here("plots/PB_dx_genes/spaghetti_plot/",
  "spatially_divergent_genes_w_all.pdf")
)
ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "SST",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "SST",
    .center = FALSE
  ),
  nrow = 2
)

ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "FKBP5",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "FKBP5",
    .center = FALSE
  ),
  nrow = 2
)

ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "SOD2",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "SOD2",
    .center = FALSE
  ),
  nrow = 2
)

ggpubr::ggarrange(
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "ALDH1A1",
    .center = TRUE
  ),
  make_spaghetti_plot(
    fc_mat = log2fc_mat,
    bg_genes = bg_gene_all,
    genes = "ALDH1A1",
    .center = FALSE
  ),
  nrow = 2
)
dev.off()


# Session Info ----
sessioninfo::session_info()