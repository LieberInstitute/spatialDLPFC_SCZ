# 80GB Memory

# Load Library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scater)
  library(tidyverse)
  # library(spatialLIBD)
  library(sessioninfo)
})

# Load Data ----
spe <- readRDS(
  here::here(
    "processed-data", "rds", "01_build_spe",
    "test_raw_spe_w_spg_N63_no_img.rds"
  )
)

# Pre-QC filtering ----
## Remove out tissue spots ----
spe <- spe[, spe$in_tissue == TRUE]
ncol(spe)
## Remove sum_umi=0 spots ----
# Remove these spots to do lognormalization
spe <- spe[, spe$sum_umi != 0]
ncol(spe)

# imgData(spe) <- NULL

# Preprocessing ----
## Fetch demo info
spe$dx <- metadata(spe)$dx_df$dx[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$sex <- metadata(spe)$dx_df$sex[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

spe$age <- metadata(spe)$dx_df$age[
  match(
    spe$sample_id,
    metadata(spe)$dx_df$sample_id
  )
]

## Log transformation -----
# Create logcounts
spe <- logNormCounts(spe)


# Sample-wise Statistics ----
## Not sure if this section is actually useful
# pdf(here("plots/02_visium_qc/pre_qc_stat.pdf"))

# plotColData(spe, "sum_umi", x = "sample_id") +
#   scale_y_log10()



# # Same plot as above as only scaling differently
# # plotColData(spe, "sum_umi", x = "sample_id") +
# #   scale_y_continuous(trans = "log2")

# plotColData(spe, "sum_gene", x = "sample_id")

# plotColData(spe, "expr_chrM_ratio", x = "sample_id") +
#   scale_y_continuous(trans = "logit")
# dev.off()

# Gene Variance Explained ----
## Note: not run well
# var_mat <- getVarianceExplained(
#   spe[,1:10000],
#   variables = c("dx", "sex", "age", "sample_id")
# )

# saveRDS(var_mat,
#         here("varExplained.rds"))



# Latent Embedding ----
## 500 HVG ----
### PCA (500 HVG) ----
spe <- runPCA(spe)

plotPCA(
  spe,
  point_size = 0.1
)

# Understand the variance explained by PCs
plotExplanatoryPCs(
  spe,
  variables = c("dx", "sex", "age", "sample_id")
)

### UMAP (500 HVG) ----

spe <- runUMAP(spe, dimred = "PCA")

library(RColorBrewer)
set.seed(1)
# Create a color palette function
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))
# Generate 65 colors
colors <- color_palette(length(unique(spe$sample_id)))[sample.int(length(unique(spe$sample_id)))]

# plotUMAP(spe,
#   point_size = 0.1,
#   colour_by = "sample_id"
# ) +
#   theme(legend.position = "none") +
#   scale_colour_manual(values = colors)

# for(.smp %in% unique(spe$sample_id)){
#   .smp <- unique(spe$sample_id)[1]
#   spe$sub_sample <- spe$sample_id == .smp

#   # Not enough memory for this shit on local computer.
#   ret_p <- plotUMAP(
#     spe[, order(spe$sub_sample)],
#    point_size = 0.1,
#    colour_by = "sub_sample"
#   ) +
#   labs(title = .smp)

#   ret_p |> print()
# }


## 2000 HVG ----
### PCA (2000 HVG) ----
set.seed(1)
spe <- runPCA(spe, ntop = 2000, name = "PCA_2000")
spe <- runUMAP(spe, dimred = "PCA_2000", name = "UMAP_2000")


## spatialDLPFC Marker Gene ----
### Find Marker Genes ----
file_DLPFC_enrich_csv <- here(
  "code/spatial_clustering/PRECAST",
  "TableS8_sig_genes_FDR5perc_enrichment.csv"
)
n_marker_gene <- 100
gene_df_raw <- read.csv(
  file_DLPFC_enrich_csv
)

gene_df <- gene_df_raw |> 
  filter(spatial_domain_resolution == "Sp09") |> 
  group_by(test) |> 
  arrange(fdr, .by_group = TRUE) |> 
  slice_head(n=n_marker_gene)


stopifnot(all(gene_df$model_type == "enrichment"))
stopifnot(nrow(gene_df) == 9*n_marker_gene)

cat("NOTE (boyiguo1): ",
    gene_df$ensembl |> unique() |> length(),
    " unique gene markers are selected for spatial clustering. \n")

# gene_df$ensembl 
# rowData(spe) |> names()
# rowData(spe)$gene_id |> head()
# sum(rowData(spe)$gene_id %in% gene_df$ensembl)

### PCA (spatialDLPFC) ----
set.seed(1)
spe <- runPCA(
  spe,
  subset_row = rowData(spe)$gene_id %in% gene_df$ensembl,
  name = "PCA_spatialDLPFC"
)

spe <- runUMAP(
  spe,
  dimred = "PCA_spatialDLPFC",
  name = "UMAP_spatialDLPFC"
)

# reducedDimNames(spe)
saveRDS(
  spe,
  here(
    "processed-data/rds/01_build_spe",
    "test_raw_spe_UMAP_N63.rds"
    )
  )

library(RColorBrewer)
set.seed(1)
# Create a color palette function
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))
# Generate 65 colors
colors <- color_palette(length(unique(spe$sample_id)))[sample.int(length(unique(spe$sample_id)))]

# Combined Plot
# png(here("plots/02_visium_qc/UMAP_2000/UMAP_2000_all_spots.png"))
# plotReducedDim(spe,
#   dimred = "UMAP_2000",
#   point_size = 0.1,
#   colour_by = "sample_id"
# ) +
#   theme(legend.position = "none") +
#   scale_colour_manual(values = colors)
# dev.off()

# Sample Specific UMAP
# for (i in 1:3){
for (i in seq.int(unique(spe$sample_id))) {
  .smp <- unique(spe$sample_id)[i]
  spe$sub_sample <- FALSE
  spe$sub_sample[spe$sample_id == .smp] <- TRUE

  # order(spe$sub_sample)

  png(
    here(
      "plots/02_visium_qc/UMAP_2000/",
      paste0("UMAP_2000_", .smp, ".png")
    )
  )
  print(
    plotReducedDim(
      spe[, order(spe$sub_sample)],
      dimred = "UMAP_2000",
      point_size = 0.1,
      colour_by = "sub_sample"
    ) +
      theme(legend.position = "none") +
      scale_colour_manual(
        values = c("FALSE" = "lightgrey", "TRUE" = colors[i])
      ) +
      labs(title = .smp)
  )
  dev.off()
}



# UMAP_2000_mat <- reducedDim(spe)
# # UMAP_2000_mat <- reducedDim(spe, "UMAP_2000")
# UMAP1_range <- c(min(UMAP_2000_mat[,1]), max(UMAP_2000_mat[,1]))
# UMAP2_range <- c(min(UMAP_2000_mat[,2]), max(UMAP_2000_mat[,2]))
# UMAP_plot_list <- list()
# pdf(here("test_umap.pdf"), height = 3, width = 3)
# for (i in 1:3) {
#   # for(i in seq.int(unique(spe$sample_id))){
#   # UMAP_plot_list <- seq.int(unique(spe$sample_id)) |>
#   # purrr::map(.f = function(i){
#   .smp <- unique(spe$sample_id)[i]
#   spe$sub_sample <- TRUE

#   # Not enough memory for this shit on local computer.
#   # ret_p <- plotUMAP(
#   print(
#     plotUMAP(
#       spe[, spe$sample_id == .smp],
#       #  point_size = 0.1,
#       colour_by = "sub_sample"
#     ) +
#       labs(title = .smp) +
#       scale_colour_manual(
#         values = c("FALSE" = "grey", "TRUE" = colors[[i]])
#       ) +
#       theme(legend.position = "none") #+
#     # scale_x_continuous(limits = UMAP1_range)+
#     # scale_y_continuous(limits = UMAP2_range)
#   )

#   # UMAP_plot_list <- c(UMAP_plot_list, ret_p)
#   # })
# }
# dev.off()


# # for (i in 1:4){
# # for(i in seq.int(unique(spe$sample_id))){
# UMAP_plot_list <- seq.int(unique(spe$sample_id)) |>
#   purrr::map(.f = function(i) {
#     .smp <- unique(spe$sample_id)[i]
#     spe$sub_sample <- spe$sample_id == .smp

#     # Not enough memory for this shit on local computer.
#     # ret_p <- plotUMAP(
#     plotUMAP(
#       spe[, spe$sample_id == .smp],
#       point_size = 0.1,
#       colour_by = "sub_sample"
#     ) +
#       labs(title = .smp) +
#       scale_colour_manual(
#         values = c("FALSE" = "grey", "TRUE" = colors[[i]])
#       ) +
#       theme(legend.position = "none") +
#       scale_x_continuous() +
#       scale_y_continuous()

#     # UMAP_plot_list <- c(UMAP_plot_list, ret_p)
#   })
# # }


# ggpubr::ggarrange(
#   plotlist = UMAP_plot_list,
#   ncol = 8
# ) |>
#   ggsave(
#     filename = here("test_UMAP_per_sample.png"),
#     plot = _
#   )

# # Understand the variance explained by PCs
# plotExplanatoryPCs(
#   spe,
#   variables = c("dx", "sex", "age", "sample_id")
# )





# plotUMAP(spe,
#   point_size = 0.1,
#   colour_by = "dx"
# ) +
#   theme(legend.position = "none")

# SpotSweeper Outliers

# outlier_spots_key <- readRDS(here("processed-data/visium_qc/outlier_key_ngb_36.rds"))

# spe$local_outliers <- "FALSE"
# spe$local_outliers[spe$key %in% outlier_spots_key] <- "TRUE"

# plotUMAP(spe,
#   point_size = 0.1,
#   colour_by = "local_outliers"
# ) +
#   theme(legend.position = "none")

# artf_spots_key <- readRDS(
#   here("processed-data/visium_qc/scratch_key.rds"),
# )
# spe$scratch <- "No"
# spe$scratch[spe$key %in% artf_spots_key] <- "Yes"


# plotUMAP(spe,
#   point_size = 0.1,
#   colour_by = "scratch"
# ) +
#   theme(legend.position = "none")






# Session ----
session_info()
