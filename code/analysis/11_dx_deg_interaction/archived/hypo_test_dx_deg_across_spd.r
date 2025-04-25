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
    "processed-data/PB_dx_genes",
    "test_inter_PRECAST_07_20240627.rds"
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

## Rownames ----
#  [1] "dxntc"              "dxscz"              "spdspd02"
#  [4] "spdspd03"           "spdspd04"           "spdspd05"
#  [7] "spdspd06"           "spdspd07"           "age"
# [10] "sexM"               "slide_idV12F14-053" "slide_idV12F14-057"
# [13] "slide_idV13F27-293" "slide_idV13F27-294" "slide_idV13F27-295"
# [16] "slide_idV13F27-296" "slide_idV13F27-336" "slide_idV13M06-279"
# [19] "slide_idV13M06-280" "slide_idV13M06-281" "slide_idV13M06-282"
# [22] "slide_idV13M06-340" "slide_idV13M06-342" "slide_idV13M06-343"
# [25] "slide_idV13M06-344" "dxscz:spdspd02"     "dxscz:spdspd03"
# [28] "dxscz:spdspd04"     "dxscz:spdspd05"     "dxscz:spdspd06"
# [31] "dxscz:spdspd07"

# Test any non-zero FC across SpD ----
contrast_fit <- contrasts.fit(fit, cont.mat)
contrast_fit <- eBayes(contrast_fit)

# Test results ----
cont_df <- topTable(contrast_fit,
  coef = sprintf("spd%02d", 1:7),
  num = Inf
) |> 
rownames_to_column(var = "ensembl")

## Create gene names ---
gene_df <- read_csv(
  here(
    "processed-data/PB_dx_genes/",
    "test_PRECAST_07.csv"
  )
) |> select(
  ensembl, gene
)

cont_df <-  cont_df |>
  left_join(
    gene_df,
    by = c("ensembl")
  )

## Save test results ----
saveRDS(
  cont_df,
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_hypo_test_non_zero_logfc_PRECAST_07.rds"
  )
)

# Estimated log2FC ----
## Format outcome
log2fc_mat <- cont_df |>
  select(
    ensembl, gene,
     starts_with("spd")
    )

saveRDS(
  log2fc_mat,
  here(
    "processed-data/PB_dx_genes/interaction",
    "test_laminar_specific_log2FC.rds"
  )
)

# Session Info ----
sessioninfo::session_info()
