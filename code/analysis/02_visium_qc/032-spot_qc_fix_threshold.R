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

raw_spe <- readRDS(
  path_raw_spe
)

spe <- raw_spe
spe <- spe[, spe$in_tissue == "TRUE"]

spe$outlier_umi <- spe$sum_umi <= 500
spe$outlier_gene <- spe$sum_gene <=200


##  Visualize Outliers ------------------------------------------------------
library(escheR) # NOTE: escheR > 1.1.1




expand.grid(
  sample_id = spe$sample_id |> unique()#,
  # method = c("scater", "OT", "joint")
) |> 
  pwalk(
    .f = function(sample_id, method){
      
      # spe_in_tissue <- spe[, spe$in_tissue == TRUE]
      
      browser()
      
      # UMI Plots
      sum_UMI_plots <- #lapply(
        # c("scater", "OT", "joint"),
        # FUN = function(method){
          make_escheR(spe[, spe$sample_id == sample_id]) |> 
            add_fill(var = "sum_umi") |> 
            add_ground(var = "outlier_umi", stroke = 0.5) +
            scale_fill_viridis_c(trans="log2") +
            scale_colour_manual(values = c("transparent", "red")) +
            labs(Title = "sum_umi")
        # }
      # ) |> 
        # ggarrange(plotlist = _, ncol = 3, nrow = 1,
        #           labels = c("scater", "OT", "joint"),
        #           common.legend = TRUE,
        #           legend = "none")
      
      sum_gene_plots <- #lapply(
        # c("scater", "OT", "joint"),
        # FUN = function(method){
          make_escheR(spe[, spe$sample_id == sample_id]) |> 
            add_fill(var = "sum_gene") |> 
            add_ground(var = "outlier_umi", stroke = 0.5) +
            scale_fill_viridis_c(trans="log2") +
            scale_colour_manual(values = c("transparent", "red")) +
            labs(Title = "sum_gene")
        # }
      # ) |> 
        # ggarrange(plotlist = _, ncol = 3, nrow = 1,
                  # labels = c("scater", "OT", "joint"),
                  # common.legend = TRUE,
                  # legend = "none")
      
      # mt_perc_plots <- lapply(
      #   c("scater", "OT", "joint"),
      #   FUN = function(method){
      #     make_escheR(spe[, spe$sample_id == sample_id]) |> 
      #       add_fill(var = "expr_chrM_ratio") |> 
      #       add_ground(var = paste0(method, "_mt_perc_outlier"), stroke = 0.5) +
      #       scale_fill_viridis_c() +
      #       scale_colour_manual(values = c("transparent", "red")) +
      #       labs(Title = "Mito %")
      #   }
      # ) |> 
      #   ggarrange(plotlist = _, ncol = 3, nrow = 1,
      #             # labels = c("scater", "OT", "joint"),
      #             common.legend = TRUE,
      #             legend = "none"
      #   )
      
      
      
      ggpubr::ggarrange(
        sum_UMI_plots,
        sum_gene_plots,
        mt_perc_plots,
        # make_escheR(spe[, spe$sample_id == sample_id]) |> 
        #   add_fill(var = "sum_umi") |> 
        #   add_ground(var = paste0(method, "_umi_outlier"), stroke = 0.5) +
        #   scale_fill_viridis_c(trans="log2") +
        #   scale_colour_manual(values = c("transparent", "red")) +
        #   labs(Title = "sum_umi"),
        # 
        # make_escheR(spe[, spe$sample_id == sample_id]) |> 
        #   add_fill(var = "sum_gene") |> 
        #   add_ground(var = paste0(method, "_gene_outlier"), stroke = 0.5) +
        #   scale_fill_viridis_c(trans="log2") +
        #   scale_colour_manual(values = c("transparent", "red")) +
        #   labs(Title = "sum_gene"),
        # 
        # make_escheR(spe[, spe$sample_id == sample_id]) |> 
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
