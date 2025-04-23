# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
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

gsea_res_df <- gsea_res_files |>
  set_names(
    # extract spd index from file names
    str_extract(
      gsea_res_files,
      "(?<=gsea_int_).*?(?=\\.csv)"
    )
  ) |>
  map(
    ~ read_csv(.x, col_types = cols())
  )


# Make upset plot ----
# TODO: move to libray session
library(UpSetR)

# Format the gsea results to acceptable format for upset plot
enrich_list <- gsea_res_df |>
  map(\(x) x$ID)

# Make upset plot
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


# save upset plot



# Session info ----
sessioninfo::session_info()
