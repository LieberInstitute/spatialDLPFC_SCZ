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
    "test_permed_laminar_specific_log2FC.rds"
  )
) |> slice_tail(n = 50)

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

  return(ret_p)
}

treg_genes <- c("AKT3", "MALAT1", "ARID1B")



make_spaghetti_plot(
  fc_mat = log2fc_mat,
  bg_genes = log2fc_mat$gene,
  genes = "",
  .center = FALSE
)

make_spaghetti_plot(
  fc_mat = log2fc_mat,
  bg_genes = log2fc_mat$gene,
  genes = "",
  .center = TRUE
)



# Session Info ----
sessioninfo::session_info()