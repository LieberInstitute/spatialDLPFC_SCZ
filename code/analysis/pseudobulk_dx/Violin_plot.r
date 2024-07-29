library(here)
library(tidyverse)
library(SingleCellExperiment)
library(sessioninfo)

# Load data ----
## Load PB data ----
pb_spe <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    "test_spe_pseudo_PRECAST_07.rds"
  )
)


## Load Genes ----
# gene_names <- readxl::read_excel(
#   here(
#     "code/analysis/pseudobulk_dx",
#     "Escher_genelist.xlsx"
#   ),
#   col_names = FALSE
# ) |> unlist()
# names(gene_names) <- NULL

gene_names <-c("ALDH1A1", "SOD2")

gene_ensembl <- rowData(pb_spe)$gene_id[match(gene_names, rowData(pb_spe)$gene_name)]


colData(pb_spe)

spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))


expression_df <- logcounts(pb_spe)[gene_ensembl, ] |>
  data.frame() |>
  rownames_to_column(var = "gene") |>
  pivot_longer(
    cols = starts_with("V"),
    names_to = c("sample_id", "spd"),
    names_pattern = "(.*)_(spd.*)"
  ) |>
  mutate(
    gene_name = rowData(pb_spe)[gene, "gene_name"],
    sample_id = str_replace(sample_id, "\\.", "-"),
    dx = metadata(pb_spe)$dx_df$dx[
      match(sample_id, metadata(pb_spe)$dx_df$sample_id)
    ]
  ) |>
  left_join(
    spd_anno_df,
    by = "spd"
  )

for (gene in gene_names) {
  # browser()
  # Create data frame
  ret_plot <- expression_df[
    which(expression_df$gene_name == gene),
  ] |>
    ggplot(aes(x = dx, y = value)) +
    geom_boxplot(aes(color = dx), show.legend = FALSE) +
    geom_jitter(size = 0.6, alpha = 0.3) +
    facet_grid(cols = vars(anno_lab)) +
    theme_light() +
    ylab("logCPM") +
    xlab("") +
    scale_color_manual(
      values = c(
        "ntc" = "blue",
        "scz" = "red"
      )
    ) +
    theme(
      strip.text.x = element_text(
        size = 12, color = "Black", face = "bold.italic"
      )
    )

  # use ggplot to make violin plot
  pdf(
    here(
      "plots/PB_dx_genes/violin_plot",
      paste0("box_plot_", gene, ".pdf")
    ),
    width = 10, height = 4
  )
  print(ret_plot)
  dev.off()
  # save pdf
}


# Session Info ----
sessioninfo::session_info()
