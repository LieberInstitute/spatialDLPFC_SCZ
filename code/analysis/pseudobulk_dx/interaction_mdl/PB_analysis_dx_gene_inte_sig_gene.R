# NOTE:
# Really need to clean up this file

# Load Libray ----
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(limma)
  # library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(pheatmap)
})

# spd_rds <- list.files(
#   here(
#     "processed-data", "rds", "layer_spd"
#   ),
#   pattern = ".rds"
# )

# Dx DE analysis ----
# pdf(
#   here(
#     "plots/PB_DE",
#     "test_dx_volcano_plot.pdf"
#   )
# )
.file <- "test_spe_pseudo_PRECAST_07.rds"
# {
# for (.file in spd_rds) {
# .file <- spd_rds[1]

.spd <- str_remove(.file, "test_spe_pseudo_") |>
  str_remove(".rds")

sce_pseudo <- readRDS(
  here(
    "processed-data", "rds", "layer_spd",
    .file
  )
)

sce_pseudo$spd <- sce_pseudo[[.spd]]


dx_mod <- model.matrix(
  ~ 0 + dx * spd + age + sex + slide_id,
  colData(sce_pseudo)
)

int_terms <- grep(":", dx_mod |> colnames(), value = TRUE)


#   dx_sig_gene <- read_csv(
#     here(
#       ,
# "test_PRECAST_07.csv"
#     )
#   ) |> filter(fdr_scz <= 0.10)



corfit <- duplicateCorrelation(
  logcounts(sce_pseudo),
  design = dx_mod,
  block = colData(sce_pseudo)$sample_id
)


# Naive without adjusting for sample-id random effect
fit <- lmFit(
  logcounts(sce_pseudo),
  design = dx_mod,
  block = colData(sce_pseudo)$sample_id,
  correlation = corfit$consensus
)
fit <- eBayes(fit)

saveRDS(
  fit,
  here(
    "processed-data/PB_dx_genes",
    "test_inter_PRECAST_07_20240627.rds"
  )
)


# Not: this is to test if interaction terms are significant, which is different to answer the question that if a there is a different for each laminar level between SCZ and NTC, or at least doesn't calculate the log fold change.
Int_df <- topTable(fit, coef = c(int_terms), num = Inf)

# cont.dif <- makeContrasts(
#   spd01 = dxscz - dxntc,
#   # spd02 = (dxscz + dxscz:spdspd02) - (dxntc),
#   levels = dx_mod
# )




cont.mat <- rbind(
  rep(-1, 7),
  rep(1, 7),
  matrix(0, nrow = 23, ncol = 7),
  cbind(rep(0, 6), diag(nrow = 6, ncol = 6))
)
#         [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#  [1,]   -1   -1   -1   -1   -1   -1   -1
#  [2,]    1    1    1    1    1    1    1
#  [3,]    0    0    0    0    0    0    0
#  [4,]    0    0    0    0    0    0    0
#  [5,]    0    0    0    0    0    0    0
#  [6,]    0    0    0    0    0    0    0
#  [7,]    0    0    0    0    0    0    0
#  [8,]    0    0    0    0    0    0    0
#  [9,]    0    0    0    0    0    0    0
# [10,]    0    0    0    0    0    0    0
# [11,]    0    0    0    0    0    0    0
# [12,]    0    0    0    0    0    0    0
# [13,]    0    0    0    0    0    0    0
# [14,]    0    0    0    0    0    0    0
# [15,]    0    0    0    0    0    0    0
# [16,]    0    0    0    0    0    0    0
# [17,]    0    0    0    0    0    0    0
# [18,]    0    0    0    0    0    0    0
# [19,]    0    0    0    0    0    0    0
# [20,]    0    0    0    0    0    0    0
# [21,]    0    0    0    0    0    0    0
# [22,]    0    0    0    0    0    0    0
# [23,]    0    0    0    0    0    0    0
# [24,]    0    0    0    0    0    0    0
# [25,]    0    0    0    0    0    0    0
# [26,]    0    1    0    0    0    0    0
# [27,]    0    0    1    0    0    0    0
# [28,]    0    0    0    1    0    0    0
# [29,]    0    0    0    0    1    0    0
# [30,]    0    0    0    0    0    1    0
# [31,]    0    0    0    0    0    0    1

# TODO: to print out
# dx_mod
# TODO: preint out contrast matrix

dx_mod |> colnames()

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

colnames(cont.mat) <- sprintf("spd%02d", 1:7)

contrast_fit <- contrasts.fit(fit, cont.mat)
contrast_fit <- eBayes(contrast_fit)






#
topTable(contrast_fit) |> View()

cont_df <- topTable(contrast_fit, coef = sprintf("spd%02d", 1:7), num = Inf)
# colnames(cont_df) <- paste0(colnames(cont_df), "_contrast")
cont_df <- cont_df |> rownames_to_column("gene_id")


gene_names_hc_ordered <- readRDS(
  here(
    "code/analysis/pseudobulk_dx",
    "spd_hierarchical_cluster_order.rds"
  )
)

all_gene_mat <- cont_df |>
  left_join(
    rowData(sce_pseudo) |> data.frame() |> select(gene_id, gene_name)
  ) |>
  filter(gene_name %in% gene_names_hc_ordered) |>
  column_to_rownames("gene_name")

stopifnot(
  nrow(all_gene_mat) == length(gene_names_hc_ordered)
)




# all_gene_mat <- cont_df |>
#   rownames_to_column(var = "gene_id") |>
#   left_join(
#     rowData(sce_pseudo) |> data.frame() |> select(gene_id, gene_name)
#   )
# rownames(all_gene_mat) <- NULL
# all_gene_mat <- all_gene_mat |> column_to_rownames("gene_name")

# all_gene_anno_row <-

spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(anno_lab = paste0(label, " (", spd, ") "))


heatmap_mat <- all_gene_mat |>
  select(starts_with("spd")) |>
  data.matrix()
colnames(heatmap_mat) <- spd_anno_df$anno_lab[match(colnames(heatmap_mat), spd_anno_df$spd)]





pdf(
  here(
    "plots/PB_dx_genes",
    "dx_sig_gene_layer_specific_heatmap.pdf"
  ),
  height = 20
)
pheatmap(
  heatmap_mat[
    gene_names_hc_ordered,
    order(colnames(heatmap_mat))
  ],
  # scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  annotation_row = all_gene_mat |> transmute(`-log10P` = -1 * log10(P.Value))
)
dev.off()


all_gene_mat |>
  filter(adj.P.Val <= 0.05) |>
  rownames()


# intersect(dx_sig_gene$gene, int_sig_df$gene_name)


# test of the deviation (FALSE)
topTable(fit,
  coef = int_terms[6] # , n=Inf
)

# Test for diagnosis
topTable(fit, coef = "dxscz", n = 30)


# TODO: add one that adjusts effect


#   dx_block_cor <-
#     registration_block_cor(
#       sce_pseudo,
#       registration_model = dx_mod,
#       var_sample_id = "sample_id"
#     )

#   dx_res <- registration_stats_enrichment(
#     sce_pseudo,
#     block_cor = dx_block_cor,
#     covars = c("spd", "age", "sex", "slide_id"),
#     var_registration = "dx",
#     gene_ensembl = "gene_id",
#     gene_name = "gene_name"
#   )

#   dir.create(
#     here("processed-data/PB_dx_genes"),
#     showWarnings = FALSE
#   )

#   ## Save DE results ----
#   dx_res |>
#     arrange(fdr_scz) |>
#     write.csv(
#       here(
#         "processed-data/PB_dx_genes",
#         paste0("test_", .spd, ".csv")
#       )
#     )

#   ## Volcano Plot ----
#   out_gene_df <- dx_res |>
#     arrange(fdr_scz) |>
#     slice_head(n = 1)
# }
# dev.off()

# Session Info ----
sessioninfo::session_info()
