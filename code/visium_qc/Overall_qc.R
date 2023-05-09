# Load Packages -----------------------------------------------------------
library(here)
library(SpatialExperiment)
library(scater)
# library(pryr)                 # Check spe size
# library(spatialLIBD)
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

dir.create(fldr_qc_plots, recursive = T)
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
### scater Threshold -----------------------------------------------------------
spe_in_tissue <- spe[, spe$in_tissue==TRUE]

scater_umi_thres <- isOutlier(spe_in_tissue$sum_umi,
                              type = "lower",
                              log = TRUE,
                              batch = spe_in_tissue$sample_id) |>
  attr(which = "thresholds")

scater_gene_thres <- isOutlier(spe_in_tissue$sum_gene,
                               type = "lower",
                               log = TRUE,
                               batch = spe_in_tissue$sample_id) |>
  attr(which = "thresholds")

scater_mt_perc_thres <- isOutlier(spe_in_tissue$expr_chrM_ratio,
                                  type = "higher",
                                  log = TRUE,
                                  batch = spe_in_tissue$sample_id) |>
  attr(which = "thresholds")

# Making sure the samples are in the same order.
stopifnot(colnames(scater_umi_thres) == colnames(scater_gene_thres),
          colnames(scater_umi_thres) == colnames(scater_mt_perc_thres))

scater_thres_df <- 
  data.frame(
    scater_umi_thres = (scater_umi_thres |> t()) [, "lower"],
    scater_gene_thres = (scater_gene_thres |> t()) [, "lower"],
    scater_mt_perc_thres = (scater_mt_perc_thres |> t()) [, "higher"]
  ) |> 
  rownames_to_column(var = "sample_id")


### Out_tissue (OT) Threshold -------------------------------------
col_df <- colData(spe) |> data.frame()

OT_thres_df <- col_df |> filter(in_tissue == FALSE) |> 
  group_by(sample_id) |> 
  summarize(
    OT_umi_thres = median(sum_umi),
    OT_gene_thres = median(sum_gene),
    OT_mt_perc_thres = median(expr_chrM_ratio)
  ) |> 
  ungroup()

### Joint (OT, scater) Threshold -------------------------------------
col_df <- col_df |> 
  left_join(scater_thres_df, by = "sample_id") |> 
  left_join(OT_thres_df, by = "sample_id") |> 
  mutate(
    # Create joint threshold, min(scater, OT)
    joint_umi_thres = min(scater_umi_thres, OT_umi_thres),
    joint_gene_thres = min(scater_gene_thres, OT_gene_thres),
    joint_mt_perc_thres = max(scater_mt_perc_thres, OT_mt_perc_thres),
    # Column for scater outliers
    scater_umi_outlier = factor(sum_umi <= scater_umi_thres),
    scater_gene_outlier = factor(sum_gene <= scater_gene_thres),
    scater_mt_perc_outlier = factor( expr_chrM_ratio >= scater_mt_perc_thres),
    # Column for OT outliers
    OT_umi_outlier = factor(sum_umi <= OT_umi_thres),
    OT_gene_outlier = factor(sum_gene <= OT_gene_thres),
    OT_mt_perc_outlier = factor( expr_chrM_ratio >= OT_mt_perc_thres),
    # Column for joint outliers
    joint_umi_outlier = factor(sum_umi <= joint_umi_thres),
    joint_gene_outlier = factor(sum_gene <= joint_gene_thres),
    joint_mt_perc_outlier = factor( expr_chrM_ratio >= joint_mt_perc_thres)
  ) # |> 
  # select(-ends_with("thres"))

# Force out_tissue to be NA
col_df[col_df$in_tissue==FALSE, endsWith(colnames(col_df), "outlier")] <- NA



colData(spe) <- DataFrame(col_df)




# Plot Thresholds ---------------------------------------------------------
# TODO: Re-write these
# metadata(spe) <- list(
#   sum_umi_thres = lib_size_thres,
#   sum_gene_thres = n_gene_thres,
#   expr_chrM_ratio_thres = mt_perc_thres
# )
# plot_metric(spe, "sum_umi")
# plot_metric(spe, "sum_gene")
# plot_metric(spe, "expr_chrM")
# plot_metric(spe, "expr_chrM_ratio", include_log = FALSE)
# # TODO:
# # plot_metric(spe, "cell_count", include_log =FALSE)

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