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

spd_lookup_tab <- data.frame(
  spd_labels = c(
    "SpD07-L1/M",
    "SpD06-L2/3",
    "SpD02-L3/4",
    "SpD05-L5",
    "SpD03-L6",
    "SpD01-WMtz",
    "SpD04-WM"
  )
) |> mutate(
  spd_raw = str_split_i(
    spd_labels, "-", 1
  ) |> tolower()
)


names(spd_deg_list) <- spd_lookup_tab$spd_labels[match(
  names(spd_deg_list), spd_lookup_tab$spd_raw
)]



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
      geom_point(size = 0.2) +
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
        x = "log2(FC)"
      ) +
      scale_color_identity(
        name = "Significance",
        labels = c(
          "FDR > 0.10",
          "Up-regulated",
          "Down-regulated"
        )
      ) +
      theme_classic(base_size = 6) +
      theme(plot.title = element_text(hjust = 0.5, size = 6, face = "bold"))
  })

## Save volcano plots ----
p_spd07_L1 <- p_volca_list[["SpD07-L1/M"]]

ggsave(
  filename = here(
    "plots/11_dx_deg_interaction",
    "volcano_plot_layer_restricted_spd07_L1_M.pdf"
  ),
  plot = p_spd07_L1,
  width = 1.13, height = 1.31,
  units = "in"
)


p_spd01 <- p_volca_list[["SpD01-WMtz"]]

ggsave(
  filename = here(
    "plots/11_dx_deg_interaction",
    "volcano_plot_layer_restricted_spd01_WMtz.pdf"
  ),
  plot = p_spd07_L1,
  width = 1.13, height = 1.31,
  units = "in"
)


# ggpubr::ggarrange(
#   plotlist = p_volca_list,
#   ncol = 5, nrow = 2
# ) |>
#   ggsave(
#     filename = here(
#       "plots/11_dx_deg_interaction",
#       "volcano_plot_layer_restricted.pdf"
#     ),
#     plot = _,
#     width = 8, height = 6
#   )


# Session Info ----
sessioninfo::session_info()
