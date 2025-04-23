# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggrepel)
})

# Load data ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

spd_deg_list <-
  spd_files |>
  set_names(
    str_extract(
      spd_files,
      "(?<=layer_specific_logFC_).*?(?=\\.csv)"
    )
  ) |>
  map(
    ~ read_csv(.x)
  )

# Make volcano plot ----
## Make a list of volcano plots ----

p_volca_list <- spd_deg_list |>
  imap(.f = function(df, spd_it) {
    # Create color categories for the volcano plot
    # browser()
    df <- df |> mutate(
      gene_cat = case_when(
        sig_10 == FALSE ~ "grey",
        sig_10 == TRUE & logFC >= 0 ~ "red",
        sig_10 == TRUE & logFC < 0 ~ "blue"
      )
    )

    ggplot(
      data = df,
      aes(x = logFC, y = -log10(P.Value), color = gene_cat)
    ) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        color = "grey"
      ) +
      geom_point() +
      # TODO:
      # Add labels on genes to higlight
      # geom_label_repel(
      # TODO: use case_when with spd_it to switch dataset.
      #   data = sig_gene_df, # Add labels last to appear as the top layer
      #   aes(label = gene),
      #   force = 2,
      #   nudge_y = 0.1
      # ) +
      labs(
        title = spd_it,
        y = "-log10(p-value)",
        x = "log2(FC) in SCZ"
      ) +
      scale_color_identity(
        name = "Significance",
        labels = c(
          "FDR > 0.10",
          "Up-regulated",
          "Down-regulated"
        )
      ) +
      theme_minimal()
  })

## Save volcano plots ----
ggpubr::ggarrange(
  plotlist = p_volca_list,
  ncol = 5, nrow = 2
) |>
  ggsave(
    filename = here(
      "plots/11_dx_deg_interaction",
      "volcano_plot_layer_specific.pdf"
    ),
    plot = _,
    width = 8, height = 6
  )


# Session Info ----
sessioninfo::session_info()
