# Load packages ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sessioninfo)
})

# Load DEG csv ----

csv_files <- list.files(
  here(
    "processed-data/PB_dx_genes"
  ),
  pattern = ".csv"
)

all_deg_df <- csv_files |>
  map_dfr(
    .f = function(.file) {
      .spd <- str_extract(.file, "PRECAST_\\d{1,2}")

      cat("Working with the result from ", .spd, "\n")
      deg_df <- read_csv(here(
        "processed-data/PB_dx_genes",
        .file
      ))

      deg_df |>
        filter(fdr_scz <= 0.05) |>
        select(ends_with("_scz"), ensembl, gene) |>
        mutate(spd = .spd)
    }
  )

# Num of DEG per SPD ----
all_deg_df |>
  group_by(spd) |>
  summarize(n = n())

# # A tibble: 14 × 2
#    spd            n
#    <chr>      <int>
#  1 PRECAST_10   248
#  2 PRECAST_11   277
#  3 PRECAST_12   178
#  4 PRECAST_13   142
#  5 PRECAST_14   241
#  6 PRECAST_15   238
#  7 PRECAST_16   194
#  8 PRECAST_3      9
#  9 PRECAST_4     34
# 10 PRECAST_5     92
# 11 PRECAST_6    121
# 12 PRECAST_7    186
# 13 PRECAST_8    223
# 14 PRECAST_9    235

# Overlapping DEGs ----
gene_sig_count <- all_deg_df |>
  group_by(gene) |>
  summarize(n = n()) |>
  arrange(desc(n))

write_csv(
  gene_sig_count,
  here(
    "processed-data/PB_dx_genes",
    "summarize_dx_gene_across_spd.csv"
  )
)


# gene_sig_porp <- gene_sig_porp |>
# mutate(n_fac = factor(n, levels = 14:1))
# table(gene_sig_porp)
# prop.table(table(gene_sig_porp$n))


# Session Info ----
session_info()
