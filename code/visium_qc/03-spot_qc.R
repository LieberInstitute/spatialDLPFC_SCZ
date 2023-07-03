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


joint_thres_df <- full_join(
  scater_thres_df, OT_thres_df, by = "sample_id"
) |> rowwise() |> 
  mutate(
    # Create joint threshold, min(scater, OT)
    joint_umi_thres = min(scater_umi_thres, OT_umi_thres),
    joint_gene_thres = min(scater_gene_thres, OT_gene_thres),
    joint_mt_perc_thres = max(scater_mt_perc_thres, OT_mt_perc_thres)
  ) |> 
  ungroup()

### Joint (OT, scater) Threshold -------------------------------------
col_df <- col_df |> 
  left_join(joint_thres_df, by = "sample_id") |> 
  mutate(
    # Column for scater outliers
    scater_umi_outlier = (sum_umi <= scater_umi_thres),
    scater_gene_outlier = (sum_gene <= scater_gene_thres),
    scater_mt_perc_outlier = ( expr_chrM_ratio >= scater_mt_perc_thres),
    # Column for OT outliers
    OT_umi_outlier = (sum_umi <= OT_umi_thres),
    OT_gene_outlier = (sum_gene <= OT_gene_thres),
    OT_mt_perc_outlier = ( expr_chrM_ratio >= OT_mt_perc_thres),
    # Column for joint outliers
    joint_umi_outlier = (sum_umi <= joint_umi_thres),
    joint_gene_outlier = (sum_gene <= joint_gene_thres),
    joint_mt_perc_outlier = ( expr_chrM_ratio >= joint_mt_perc_thres)
  ) # |> 
# select(-ends_with("thres"))

# Force out_tissue to be NA
col_df[col_df$in_tissue==FALSE, endsWith(colnames(col_df), "outlier")] <- NA


# Summary Statistics of Outlier Spots -------------------------------------
## scater -------------------------------------
col_df |> 
  filter(in_tissue == TRUE) |> 
  group_by(sample_id) |> 
  summarize(
    n_umi_outlier = sum(scater_umi_outlier),
    perc_umi_outlier = n_umi_outlier/n(),
    n_gene_outlier = sum(scater_umi_outlier),
    perc_gene_outlier = n_gene_outlier/n(),
    n_mt_perc_outlier = sum(scater_mt_perc_outlier),
    perc_mt_perc_outlier = n_mt_perc_outlier/n()
  )





## OT -------------------------------------
col_df |> 
  filter(in_tissue == TRUE) |> 
  group_by(sample_id) |> 
  summarize(
    n_umi_outlier = sum(OT_umi_outlier),
    perc_umi_outlier = n_umi_outlier/n(),
    n_gene_outlier = sum(OT_umi_outlier),
    perc_gene_outlier = n_gene_outlier/n(),
    n_mt_perc_outlier = sum(OT_mt_perc_outlier),
    perc_mt_perc_outlier = n_mt_perc_outlier/n()
  )


## Joint -------------------------------------
col_df |> 
  filter(in_tissue == TRUE) |> 
  group_by(sample_id) |> 
  summarize(
    n_umi_outlier = sum(joint_umi_outlier),
    perc_umi_outlier = n_umi_outlier/n(),
    n_gene_outlier = sum(joint_umi_outlier),
    perc_gene_outlier = n_gene_outlier/n(),
    n_mt_perc_outlier = sum(joint_mt_perc_outlier),
    perc_mt_perc_outlier = n_mt_perc_outlier/n()
  )

# spe_in_tissue$low_lib_size <- qcfilter$low_lib_size |> factor()
# spe_in_tissue$low_n_features <- qcfilter$low_n_features |> factor()
# spe_in_tissue$high_subsets_Mito_percent <- qcfilter$high_subsets_Mito_percent |>
#   factor()
# Update Spe with
colData(spe) <- DataFrame(col_df)







## Plot Thresholds ---------------------------------------------------------
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
# spe[,spe$expr_chrM==0]



##  Visualize Outliers ------------------------------------------------------

# qcfilter <- data.frame(
#   low_lib_size = isOutlier(spe_in_tissue$sum_umi, type = "lower",
#                            log = TRUE, batch = spe_in_tissue$sample_id),
#   low_n_features = isOutlier(spe_in_tissue$sum_gene, type = "lower", log = TRUE,
#                              batch = spe_in_tissue$sample_id),
#   high_subsets_Mito_percent = isOutlier(spe_in_tissue$expr_chrM_ratio,
#                                         log = FALSE,
#                                         type = "higher",
#                                         batch = spe_in_tissue$sample_id)
# )


library(escheR) # NOTE: escheR > 1.1.1
library(ggpubr)
# smp_id <- "Br5367_D1"



# Depreciate
# file.path(fldr_outlier_plots,
#           c("scater", "OT", "joint")) |> 
#   walk( .f = dir.create)


#TODO: may be it would be better just make 3*4 panels
expand.grid(
  sample_id = col_df$sample_id |> unique()#,
  # method = c("scater", "OT", "joint")
) |> 
  pwalk(
    .f = function(sample_id, method){
      
      
      spe_in_tissue <- spe[, spe$in_tissue == TRUE]
      # method <- "OT"
      method = c("scater", "OT", "joint")
      
      # browser()
      
      # UMI Plots
      sum_UMI_plots <- lapply(
        c("scater", "OT", "joint"),
        FUN = function(method){
          make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == sample_id]) |> 
            add_fill(var = "sum_umi") |> 
            add_ground(var = paste0(method, "_umi_outlier"), stroke = 0.5) +
            scale_fill_viridis_c(trans="log2") +
            scale_colour_manual(values = c("transparent", "red")) +
            labs(Title = "sum_umi")
        }
      ) |> 
        ggarrange(plotlist = _, ncol = 3, nrow = 1,
                  labels = c("scater", "OT", "joint"),
                  common.legend = TRUE,
                  legend = "none")
      
      sum_gene_plots <- lapply(
        c("scater", "OT", "joint"),
        FUN = function(method){
          make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == sample_id]) |> 
            add_fill(var = "sum_gene") |> 
            add_ground(var = paste0(method, "_gene_outlier"), stroke = 0.5) +
            scale_fill_viridis_c(trans="log2") +
            scale_colour_manual(values = c("transparent", "red")) +
            labs(Title = "sum_gene")
        }
      ) |> 
        ggarrange(plotlist = _, ncol = 3, nrow = 1,
                  # labels = c("scater", "OT", "joint"),
                  common.legend = TRUE,
                  legend = "none")
      
      mt_perc_plots <- lapply(
        c("scater", "OT", "joint"),
        FUN = function(method){
          make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == sample_id]) |> 
            add_fill(var = "expr_chrM_ratio") |> 
            add_ground(var = paste0(method, "_mt_perc_outlier"), stroke = 0.5) +
            scale_fill_viridis_c() +
            scale_colour_manual(values = c("transparent", "red")) +
            labs(Title = "Mito %")
        }
      ) |> 
        ggarrange(plotlist = _, ncol = 3, nrow = 1,
                  # labels = c("scater", "OT", "joint"),
                  common.legend = TRUE,
                  legend = "none"
                  )
      
      
      
      ggpubr::ggarrange(
        sum_UMI_plots,
        sum_gene_plots,
        mt_perc_plots,
        # make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == sample_id]) |> 
        #   add_fill(var = "sum_umi") |> 
        #   add_ground(var = paste0(method, "_umi_outlier"), stroke = 0.5) +
        #   scale_fill_viridis_c(trans="log2") +
        #   scale_colour_manual(values = c("transparent", "red")) +
        #   labs(Title = "sum_umi"),
        # 
        # make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == sample_id]) |> 
        #   add_fill(var = "sum_gene") |> 
        #   add_ground(var = paste0(method, "_gene_outlier"), stroke = 0.5) +
        #   scale_fill_viridis_c(trans="log2") +
        #   scale_colour_manual(values = c("transparent", "red")) +
        #   labs(Title = "sum_gene"),
        # 
        # make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == sample_id]) |> 
        #   add_fill(var = "expr_chrM_ratio") |> 
        #   add_ground(var = paste0(method, "_mt_perc_outlier"), stroke = 0.5) +
        #   scale_fill_viridis_c() +
        #   scale_colour_manual(values = c("transparent", "red")) +
        #   labs(Title = "Mito %"),
        ncol = 1,
        labels = c("UMI", "GENE", "Mt_Perc")
      ) |> 
        ggsave(
          filename = file.path(
            fldr_outlier_plots,
            paste0(sample_id,".pdf")
          ),
          height = 12,
          width = 12
        )
    }
  )


# Remove Outliers ---------------------------------------------------------

spe$is_outlier <- spe$joint_umi_outlier | spe$joint_gene_outlier

ret_spe <- spe[, spe$in_tissue == TRUE & (!spe$is_outlier)]

# Validation
stopifnot(all(ret_spe$in_tissue ==TRUE))
stopifnot(all(ret_spe$joint_umi_outlier ==FALSE))
stopifnot(all(ret_spe$joint_gene_outlier ==FALSE))

# Save output


saveRDS(
  ret_spe,
  here::here("processed-data", "rds", "spe", "spe_after_spot_qc.rds")
)

# TODO: is there a big difference between log2 and log10 outlier detection?


# Outlier Statistics ------------------------------------------------------




# Remove Outlier Spots ----------------------------------------------------
n_outlier_spot <- sum(spe$is_outlier & spe$in_tissue, na.rm = TRUE)
n_total_spot <- sum(spe$in_tissue, na.rm = TRUE)
per_outlier_spot <- n_outlier_spot/n_total_spot

col_df <- colData(spe) |> data.frame()
per_sample_ouliter_n <- col_df |> group_by(sample_id) |> 
  summarize(n_ouliter = sum(is_outlier & in_tissue, na.rm = TRUE))

per_sample_ouliter_n
#   sample_id     n_ouliter
#    <glue>            <int>
# 1 V12F14-053_A1        53
# 2 V12F14-053_B1        99
# 3 V12F14-053_C1       112
# 4 V12F14-053_D1        74
# 5 V12F14-057_A1        71
# 6 V12F14-057_B1        31
# 7 V12F14-057_C1       101
# 8 V12F14-057_D1        74

# Median
median(per_sample_ouliter_n$n_ouliter)




# * Remove spots without counts ---------------------------------------------
# no_expr_spot <- colSums(counts(raw_spe)) != 0
# 
# spe <- spe[, no_expr_spot]
# dim(spe)


# * Remove spots outside of tissue area ------------------------------------
# stopifnot(is.logical(spe$in_tissue))
# 
# # Evaluate if in-tissue definition is accurate
# vis_grid_clus(
#   spe = spe,
#   clustervar = "in_tissue",
#   pdf_file = file.path(
#     fldr_qc_plots,
#     "in_tissue_grid_raw.pdf"  # Plot File Name
#   ),
#   sort_clust = FALSE,
#   colors = c("FALSE" = "grey90", "TRUE" = "orange")
# )
# 
# # Note, if in_tissue is integer, please convert it to a logic variable
# spe <- spe[, spe$in_tissue]
# ncol(spe)
# 
# # Visual Validation
# vis_grid_clus(
#   spe = spe,
#   clustervar = "in_tissue",
#   pdf_file = file.path(
#     fldr_qc_plots,
#     "in_tissue_grid.pdf"  # Plot File Name
#   ),
#   sort_clust = FALSE,
#   colors = c("FALSE" = "grey90", "TRUE" = "orange")
# )