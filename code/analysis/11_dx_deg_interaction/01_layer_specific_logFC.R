# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ComplexHeatmap)
  library(SingleCellExperiment)
  library(limma)
  library(sessioninfo)
  library(ggrepel)
})

# Load interaction fit ----
fit <- readRDS(
  here(
    "processed-data/rds/11_dx_deg_interaction",
    "limma_obj_int_PRECAST_07.rds"
  )
)


# Create contrast matrix of SpDs ----
cont.mat <- rbind(
  rep(-1, 7),
  rep(1, 7),
  matrix(0, nrow = 23, ncol = 7),
  cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
)
colnames(cont.mat) <- sprintf("spd%02d", 1:7)

## Output of contrasts ----
#       spd01 spd02 spd03 spd04 spd05 spd06 spd07
#  [1,]    -1    -1    -1    -1    -1    -1    -1
#  [2,]     1     1     1     1     1     1     1
#  [3,]     0     0     0     0     0     0     0
#  [4,]     0     0     0     0     0     0     0
#  [5,]     0     0     0     0     0     0     0
#  [6,]     0     0     0     0     0     0     0
#  [7,]     0     0     0     0     0     0     0
#  [8,]     0     0     0     0     0     0     0
#  [9,]     0     0     0     0     0     0     0
# [10,]     0     0     0     0     0     0     0
# [11,]     0     0     0     0     0     0     0
# [12,]     0     0     0     0     0     0     0
# [13,]     0     0     0     0     0     0     0
# [14,]     0     0     0     0     0     0     0
# [15,]     0     0     0     0     0     0     0
# [16,]     0     0     0     0     0     0     0
# [17,]     0     0     0     0     0     0     0
# [18,]     0     0     0     0     0     0     0
# [19,]     0     0     0     0     0     0     0
# [20,]     0     0     0     0     0     0     0
# [21,]     0     0     0     0     0     0     0
# [22,]     0     0     0     0     0     0     0
# [23,]     0     0     0     0     0     0     0
# [24,]     0     0     0     0     0     0     0
# [25,]     0     0     0     0     0     0     0
# [26,]     0     1     0     0     0     0     0
# [27,]     0     0     1     0     0     0     0
# [28,]     0     0     0     1     0     0     0
# [29,]     0     0     0     0     1     0     0
# [30,]     0     0     0     0     0     1     0
# [31,]     0     0     0     0     0     0     1

# Error prevention
# Design matrix is conformable with the contrast matrix
stopifnot(nrow(cont.mat) == ncol(fit$design))

# Layer specific logFC and hypo test ----

colnames(cont.mat) |>
  set_names() |>
  iwalk(.f = function(spd_it, idx) {
    # browser()
    spd_cont_df <- topTable(
      contrast_fit,
      coef = spd_it, num = Inf
    ) |>
      rownames_to_column("gene_id") |>
      mutate(
        sig_10 = adj.P.Val <= 0.1
      )

    write_csv(
      spd_cont_df,
      here(
        "processed-data/rds/11_dx_deg_interaction",
        paste0(
          "layer_specific_logFC_",
          gsub("[^[:alnum:]_]", "_", idx),
          ".csv"
        )
      )
    )
  })





