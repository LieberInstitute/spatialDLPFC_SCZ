# NOTE (boyi):
# I think this analysis would be less informative comparing to enrichment test
# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  # library(spatialLIBD)
  library(readxl)
  library(sessioninfo)
})


# Load data ----
## Load Ruzicka cell type DEGs ----
# List all sheets in the Ruzicka Excel file
ruzicka_sheets <- excel_sheets(
  here(
    "code/analysis/12_deg_integration",
    "ruzicka_cell_type_DEG.xlsx"
  )
)

# Error prevention
stopifnot(
  length(ruzicka_sheets) == 25
)

# format to list of data frames for each cell types
ruzicka_deg_list <- ruzicka_sheets |>
  set_names() |>
  map(
    ~ read_xlsx(
      here(
        "code/analysis/12_deg_integration",
        "ruzicka_cell_type_DEG.xlsx"
      ),
      sheet = .x
    ) |>
      dplyr::select(gene, starts_with("Meta_")) |>
      mutate(
        cell_type = .x,
        ruzicka_sig_gene = case_when(
          Meta_adj.P.Val < 0.05 & abs(Meta_logFC) > 0.1 ~ TRUE,
          TRUE ~ FALSE
        )
      )
  )

# Some statistics
## N genes per cell type
ruzicka_deg_list |>
  map_int(
    ~ nrow(.x)
  )

## Numb of significant genes per cell type
ruzicka_deg_list |>
  map_int(
    ~ sum(.x$ruzicka_sig_gene, na.rm = TRUE)
  )


ruzicka_total_sig_gene_n <- ruzicka_deg_list |>
  map_int(
    ~ sum(.x$ruzicka_sig_gene, na.rm = TRUE)
  ) |>
  sum()

# Error prevention
stopifnot(
  ruzicka_total_sig_gene_n == 6634
)


## Load Layer specific ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  # pattern = "layer_specific_logFC_.*\\.csv",
  pattern = "layer_restricted_logFC_.*\\.csv",
  full.names = TRUE
)

names(spd_files) <- str_extract(
  spd_files,
  # "(?<=layer_specific_logFC_).*?(?=\\.csv)"
  "(?<=layer_restricted_logFC_).*?(?=\\.csv)"
)

spd_deg_list <- map(
  spd_files,
  ~ read_csv(.x, col_types = cols()) |>
    mutate(gene = AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ))
)


# Create merged data ----
ruzick_merged_df <- ruzicka_deg_list |>
  map(
    ~ .x |>
      inner_join(
        # TODO: for loop for this step, iterating over spd_deg_list
        spd_deg_list[[1]],
        by = "gene"
      )
  )


## Quick Statsitics ----
# Number of overlapping genes
ruzick_merged_df |>
  map_int(
    ~ nrow(.x)
  )

# Correlation with cell types (t_stat)
ruzick_merged_df |>
  map_dbl(
    ~ cor(
      .x$Meta_tstat,
      .x$t,
      use = "pairwise.complete.obs"
    )
  )

# Correlation with cell types (logFC)
ruzick_merged_df |>
  map_dbl(
    ~ cor(
      .x$Meta_logFC,
      .x$logFC_scz,
      use = "pairwise.complete.obs"
    )
  )

# Scatter plot ----
## Create the plot
scatter_plot <- ruzick_merged_df |>
  imap(
    ~ ggplot(
      # TODO: move this step to the data merging step
      .x |> mutate(
        study_sig = case_when(
          # TODO: confirm how is the study defines the significance
          Meta_adj.P.Val < 0.10 & fdr_scz < 0.10 ~ "Both",
          Meta_adj.P.Val >= 0.10 & fdr_scz < 0.10 ~ "LIBD_PNN",
          Meta_adj.P.Val < 0.10 & fdr_scz >= 0.10 ~ "Ruzicka",
          TRUE ~ "Neither"
        )
      ),
      aes(x = Meta_logFC, y = logFC_scz, color = study_sig)
    ) +
      geom_point(size = 0.5) +
      labs(
        title = .y,
        x = "Ruzicka logFC",
        y = "SCZ logFC"
      ) +
      theme_minimal()
  )

## Format the plot ----
scatter_plot <- scatter_plot |>
  map(~ .x +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(
      values = c("Both" = "black", "LIBD_PNN" = "blue", "Ruzicka" = "red", "Neither" = "grey")
    ))

## Save the plot ----
pdf(
  here(
    "code/analysis/12_deg_integration",
    "test_scatter_ruzick_172.pdf"
  )
)
walk(
  scatter_plot,
  ~ print(.x)
)
dev.off()
