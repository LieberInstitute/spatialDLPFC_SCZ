# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(UpSetR)
  library(sessioninfo)
})

# Load data ----
# Load all gsea results
gsea_res_files <- list.files(
  here("processed-data/rds/11_dx_deg_interaction"),
  pattern = "gsea_int_.*\\.csv",
  full.names = TRUE
)

# error prevention
stopifnot(gsea_res_files |> length() == 7)

# Load data ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

## Load spd annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(
    anno_lab = factor(
      paste0(gsub("spd", "SpD", spd), "-", label),
      levels = c(
        "SpD07-L1/M",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      )
    )
  )



gsea_res_df <- gsea_res_files |>
  set_names(
    nm = spd_anno_df[
      match(
        # extract spd index from file names
        str_extract(
          gsea_res_files,
          "(?<=gsea_int_).*?(?=\\.csv)"
        ),
        spd_anno_df$spd
      ),
      "anno_lab",
      drop = TRUE
    ]
  ) |>
  map(
    ~ read_csv(.x, col_types = cols())
  )


# Make upset plot ----
# Format the gsea results to acceptable format for upset plot
enrich_list <- gsea_res_df |>
  map(\(x) x$ID)

# Make upset plot
pdf(
  here(
    "plots/11_dx_deg_interaction",
    "upset_layer_GO_GSEA.pdf"
  ),
  width = 8,
  height = 6
)
upset(
  fromList(enrich_list),
  nsets = 7,
  # sets = sort(
  #   names(enrich_list),
  #   decreasing = TRUE
  # ),
  keep.order = TRUE,
  nintersects = NA,
  order.by = c("freq"),
  # main.bar.color = "blue",
  # sets.bar.color = "red",
  # mainbar.y.label = "Enriched Genes",
  sets.x.label = "Spatial Domains"
)
dev.off()

# Idenitfy Layer specific GO ----
unique_GO <- map(
  names(enrich_list),
  ~ setdiff(
    enrich_list[[.x]],
    Reduce(union, enrich_list[names(enrich_list) != .x])
  )
) |>
  set_names(names(enrich_list))

unique_genes_fnl <- imap_dfr(
  unique_GO,
  .f = function(.x, idx) {
    # browser()
    gsea_res_df[[idx]] |>
      dplyr::filter(ID %in% .x) |>
      select(
        ID, Description, NES,
        core_enrichment_symbol
      ) |>
      mutate(
        spd = idx
      )
  }
)

write_csv(
  unique_genes_fnl,
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "unique_layer_GO_GSEA.csv"
  )
)

# Session info ----
sessioninfo::session_info()
