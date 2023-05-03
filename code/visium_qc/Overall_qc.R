# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(scater)
# library(pryr)                 # Check spe size
library(spatialLIBD)
library(tidyverse)



# File Paths --------------------------------------------------------------
path_raw_spe <- here("processed-data/rds/spe",
                     "01_build_spe/", "spe_raw.rds")

path_clean_spe <- here("process-data/rds/spe",
                       "spe_clean.rds")

fldr_qc_plots <- here("plots", "02_visium_qc")

fldr_tissue_plots <- here("plots", "02_visium_qc",
                          "in_tissue_plots")
fldr_outlier_plots <- here("plots", "02_visium_qc",
                           "outlier_highlight")

dir.create( fldr_qc_plots, recursive = T)
dir.create(fldr_tissue_plots, recursive = T)
dir.create(fldr_outlier_plots, recursive = T)


# Source Helper Functions -------------------------------------------------

source(here("code/visium_qc/fun_plot_metrics.R"))


# QC Code -----------------------------------------------------------------
raw_spe <- readRDS(
  path_raw_spe
)

spe <- raw_spe

# is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
# # This is eqquivalent to 
# # is_mito <- which(seqnames(spe) == "chrM")
# rowData(spe)$gene_name[is_mito]           # Show mt gene names
# 
# spe <- addPerCellQC(spe, subsets = list(mito = is_mito))


## In Tissue Plot -----------------------------------------------------------------
spe$sample_id |> unique() |> 
  walk(
    .f = ~vis_clus(
      spe,
      sampleid = .x,
      clustervar = "in_tissue",
      colors = c(
        "TRUE" = "transparent",
        "FALSE" = "grey90"#, "TRUE" = "orange"
        # TODO: try this with transparent
      ) #,
      # alpha = 0.5
    ) |> 
      ggsave(
        filename = file.path(
          fldr_tissue_plots,
          paste0(.x,".pdf")
        ),
        height = 8,
        width = 8
      )
  )

## Create Thresholds -----------------------------------------------------------
spe_in_tissue <- spe[, spe$in_tissue==TRUE]

lib_size_thres <- isOutlier(spe_in_tissue$sum_umi,
                            type = "lower",
                            log = TRUE,
                            batch = spe_in_tissue$sample_id) |> 
  attr(which = "thresholds")

n_gene_thres <- isOutlier(spe[, spe$in_tissue==TRUE]$sum_gene,
                          type = "lower",
                          log = TRUE,
                          batch = spe_in_tissue$sample_id) |> 
  attr(which = "thresholds")

mt_perc_thres <- isOutlier(spe[, spe$in_tissue==TRUE]$expr_chrM_ratio,
                           type = "higher",
                           log = TRUE,
                           batch = spe_in_tissue$sample_id) |> 
  attr(which = "thresholds")


metadata(spe) <- list(
  sum_umi_thres = lib_size_thres,
  sum_gene_thres = n_gene_thres,
  expr_chrM_ratio_thres = mt_perc_thres
)


plot_metric(spe, "sum_umi")
plot_metric(spe, "sum_gene")
plot_metric(spe, "expr_chrM")
plot_metric(spe, "expr_chrM_ratio", include_log = FALSE)
# TODO:
# plot_metric(spe, "cell_count", include_log =FALSE)

# One spot Mt counts are just 0
sum(spe$expr_chrM==0)
spe[,spe$expr_chrM==0]


##  Visualize Outliers ------------------------------------------------------

qcfilter <- data.frame(
  low_lib_size = isOutlier(spe_in_tissue$sum_umi, type = "lower",
                           log = TRUE, batch = spe_in_tissue$sample_id),
  low_n_features = isOutlier(spe_in_tissue$sum_gene, type = "lower", log = TRUE,
                             batch = spe_in_tissue$sample_id),
  high_subsets_Mito_percent = isOutlier(spe_in_tissue$expr_chrM_ratio,
                                        log = FALSE,
                                        type = "higher",
                                        batch = spe_in_tissue$sample_id)
)

## Summary statistics
# Number of spots excluding
cbind(qcfilter, sample_id = spe_in_tissue$sample_id) |>
  group_by(sample_id) |> 
  summarize(
    n_lib_size = sum(low_lib_size),
    perc_lib_size = n_lib_size/n(),
    n_low_n_features = sum(low_n_features),
    perc_low_n_features = n_low_n_features/n(),
    high_subsets_Mito_percent = sum(high_subsets_Mito_percent),
  )


library(escheR) # NOTE: escheR > 0.99.8

# smp_id <- "Br5367_D1"

spe_in_tissue$low_lib_size <- qcfilter$low_lib_size |> factor()
spe_in_tissue$low_n_features <- qcfilter$low_n_features |> factor()
spe_in_tissue$high_subsets_Mito_percent <- qcfilter$high_subsets_Mito_percent |>
  factor()


spe$sample_id |> unique() |> 
  walk(
    .f = function(smp_id){
      
      ggpubr::ggarrange(
        
        make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == smp_id]) |> 
          add_fill(var = "sum_umi") |> 
          add_ground(var = "low_lib_size", stroke = 0.5) +
          scale_fill_viridis_c(trans="log2") +
          scale_colour_manual(values = c("transparent", "red")) +
          labs(Title = "Library Size"),
        
        make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == smp_id]) |> 
          add_fill(var = "sum_gene") |> 
          add_ground(var = "low_n_features", stroke = 0.5) +
          scale_fill_viridis_c(trans="log2") +
          scale_colour_manual(values = c("transparent", "red")) +
          labs(Title = "# of Features"),
        
        make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == smp_id]) |> 
          add_fill(var = "expr_chrM_ratio") |> 
          add_ground(var = "high_subsets_Mito_percent", stroke = 0.5) +
          scale_fill_viridis_c() +
          scale_colour_manual(values = c("transparent", "red")) +
          labs(Title = "Mito %"),
        ncol = 1
      ) |> 
        ggsave(
          filename = file.path(
            fldr_outlier_plots,
            paste0(smp_id,".pdf")
          ),
          height = 12,
          width = 6
        )
      
      
    }
  )

# TODO: is there a big difference between log2 and log10 outlier detection?














# * Remove spots without counts ---------------------------------------------
no_expr_spot <- colSums(counts(raw_spe)) != 0

spe <- spe[, no_expr_spot]
dim(spe)


# * Remove spots outside of tissue area ------------------------------------
stopifnot(is.logical(spe$in_tissue))

# Evaluate if in-tissue definition is accurate
vis_grid_clus(
  spe = spe,
  clustervar = "in_tissue",
  pdf_file = file.path(
    fldr_qc_plots,
    "in_tissue_grid_raw.pdf"  # Plot File Name
  ),
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "orange")
)

# Note, if in_tissue is integer, please convert it to a logic variable
spe <- spe[, spe$in_tissue]
ncol(spe)

# Visual Validation
vis_grid_clus(
  spe = spe,
  clustervar = "in_tissue",
  pdf_file = file.path(
    fldr_qc_plots,
    "in_tissue_grid.pdf"  # Plot File Name
  ),
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "orange")
)





# Feature/Gene Analysis -----------------------------------------------------
# * Remove Genes with 0 count ---------------------------------------------

gene_sum <- rowSums(counts(raw_spe))
mean(gene_sum == 0)     # High percentage of non expressing gene
spe <- raw_spe[gene_sum != 0,]
dim(spe)

# * (Optional) Remove Genes with 0 variance ---------------------------------------------
gene_var <- counts(spe) |> MatrixGenerics::rowVars()
if(mean(gene_var == 0) != 0){
  spe <- spe[gene_var != 0, ]
}


# summary(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0    66.0   170.0   291.6   326.0 39053.0

# hist(log10(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)]))
# 0  1  2  3  4  5
# 12 17 50 59 65 81


# TODO: error
# vis_grid_gene(
#   spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
#   geneid = "sum_umi",
#   pdf = here::here("plots", "01_build_spe", "out_tissue_sum_umi_all.pdf"),
#   assayname = "counts",
#   auto_crop = FALSE
# )
# # Error in frame_lims$y_min:frame_lims$y_max : NA/NaN argument
# # Auto_crop doesn't work in this case.
# 
# summary(spe_raw$expr_chrM_ratio[which(!colData(spe_raw)$in_tissue)])
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# #  0.0000  0.1500  0.1947  0.2072  0.2500  1.0000      12
# 
# vis_grid_gene(
#   spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
#   geneid = "expr_chrM_ratio",
#   pdf = here::here("plots", "01_build_spe", "out_tissue_expr_chrM_ratio_all.pdf"),
#   assayname = "counts",
#   auto_crop = FALSE
# )
# 
# vis_grid_gene(
#   spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
#   geneid = "sum_gene",
#   pdf = here::here("plots", "01_build_spe", "out_tissue_sum_gene_all.pdf"),
#   assayname = "counts"
# )
# 
# summary(spe$sum_umi)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 1    1358    2310    2674    3540   46881
# vis_grid_gene(
#   spe = spe,
#   geneid = "sum_umi",
#   pdf = here::here("plots", "01_build_spe", "in_tissue_sum_umi_all.pdf"),
#   assayname = "counts"
# )
# 
# summary(spe$sum_gene)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 1     910    1426    1508    1999    8344
# 
# vis_grid_gene(
#   spe = spe,
#   geneid = "sum_gene",
#   pdf = here::here("plots", "01_build_spe", "in_tissue_sum_gene_all.pdf"),
#   assayname = "counts"
# )
# 
# summary(spe$expr_chrM_ratio)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 0.00000 0.07279 0.10440 0.12106 0.15513 1.00000
# vis_grid_gene(
#   spe = spe,
#   geneid = "expr_chrM_ratio",
#   pdf = here::here("plots", "01_build_spe", "in_tissue_expr_chrM_ratio_all.pdf"),
#   assayname = "counts"
# )
# 
# vis_grid_gene(
#   spe = spe_raw,
#   geneid = "sum_umi",
#   pdf = here::here("plots", "01_build_spe", "all_sum_umi.pdf"),
#   assayname = "counts"
# )
# 
# vis_grid_gene(
#   spe = spe_raw,
#   geneid = "sum_gene",
#   pdf = here::here("plots", "01_build_spe", "all_sum_gene.pdf"),
#   assayname = "counts"
# )
# 
# vis_grid_gene(
#   spe = spe_raw,
#   geneid = "expr_chrM_ratio",
#   pdf = here::here("plots", "01_build_spe", "all_expr_chrM_ratio.pdf"),
#   assayname = "counts"
# )
# 
# 
# 
# # Genes with no reads -----------------------------------------------------
# stopifnot()
# 
# 
# 
# pryr::object_size(spe)

