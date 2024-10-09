# Load packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(glue)
  library(sessioninfo)
})

# Loop over all SPG + analysis level ----
# p_cut_off <- 0.10
analysis_type <- c("SPD", "donor")
spg_names <- c("pnn_pos", "neuropil_pos", "neun_pos", "vasc_pos")

analysis_combo <- expand.grid(
  type = analysis_type,
  spg = spg_names
)

pb_files <- analysis_combo |> glue_data("test_{type}_pseudo_{spg}.csv")

pb_files |>
  walk(.f = function(.pb_file) {
    # Error prevention
    # Check if pseudobulk rds file exists
    stopifnot(
      file.exists(
        here(
          "processed-data/spg_pb_de",
          .pb_file
        )
      )
    )

    ## Load DE CSV ----
    gene_df <- read_csv(
      here(
        "processed-data/spg_pb_de",
        .pb_file |> str_replace(".rds", ".csv")
      )
    ) |> mutate(
      sig_05 = fdr_scz <= 0.05,
      sig_10 = fdr_scz <= 0.10
    )



    ## Create Volcano Plot ----
    pdf(
      here(
        "plots/dx_deg_spg_pos",
        .pb_file |> str_replace(".csv", ".pdf")
      )
    )
    
    ret_p <- ggplot(
      data = gene_df,
      aes(x = logFC_scz, y = -log10(p_value_scz), color = sig_10)
    ) +
      geom_point() +
      ### (Option) Annotate Genes if needed ----
      # geom_label_repel(
      #   data = sig_gene_df, # Add labels last to appear as the top layer
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
        title = .pb_file,
        y = "-log10(p-value)",
        x = "log2(FC)"
      ) +
      scale_color_manual(
        values = c(
          "FALSE" = "grey",
          "TRUE" = "red"
        ),
        breaks = c("TRUE", "FALSE"),
        name = "FDR < 0.1"
      ) +
      theme_minimal()

    ret_p |> print()
    dev.off()
  })

# Session Info ----
session_info()
