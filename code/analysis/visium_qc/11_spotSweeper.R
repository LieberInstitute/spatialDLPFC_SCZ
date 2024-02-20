
# Load Packages -----------------------------------------------------------
library(SpatialExperiment)
library(SpotSweeper)
library(sessioninfo)
library(escheR)
library(here)
# library(scater)
# library(spatialLIBD)
# library(tidyverse)

# File Paths --------------------------------------------------------------
path_raw_spe <- here(
  "processed-data/rds",
  "01_build_spe",
  "raw_spe_wo_SPG_N63.rds")

# QC Code -----------------------------------------------------------------
raw_spe <- readRDS(
  path_raw_spe
)

spe <- raw_spe




# Find Tissue Artifacts (Tears and Hangneils)---------------------------------------------------

sample_info <- readxl::read_xlsx(
  here("raw-data/experiment_info",
       "PNN_64_collection_Summary_final.xlsx")
)

artf_samples <- sample_info |> 
  select(SlideSerial, CaptureArea, hangnail, tear) |> 
  filter(hangnail=="yes" | tear=="yes") |> 
  transmute(sample_id = paste0(SlideSerial, "_", CaptureArea)) |>
  pull(sample_id)


pdf(
  here("plots/02_visium_qc/spotSweeper", "artifacts_hangneil+tears.pdf"),
    width = 5.9, height = 3.9
  )
for(.sample in artf_samples){
  # Samples with tears
  # .sample <- artf_samples[3]
  # .sample <- "V13M06-282_B1"
  sub_spe <- spe[, spe$sample_id==.sample]
  
  sub_spe <- sub_spe[, sub_spe$in_tissue == 1]
  # features <- c('sum_umi' ,'sum_gene', 'm')
  # sub_spe<- localOutliers(
  #   sub_spe, 
  #   features=features,
  #   n_neighbors=18, 
  #   data_output=TRUE,
  #   method="multivariate"
  # )
  # 
  # make_escheR(sub_spe) |> 
  #   add_fill(var = "sum_umi_log2", point_size=1.25) |>
  #   add_ground(var = "local_outliers",point_size=1.25, stroke = 1) +
  #   scale_color_manual(
  #     name = "", # turn off legend name for ground_truth
  #     values = c(
  #       "TRUE" = "red",
  #       "FALSE" = "transparent")
  #   ) +
  #   scale_fill_gradient(low ="white",high =  "darkgreen")
  # 
  # make_escheR(sub_spe) |> 
  #   add_fill(var = "sum_gene_log2", point_size=1.25) |>
  #   add_ground(var = "local_outliers",point_size=1.25, stroke = 1) +
  #   scale_color_manual(
  #     name = "", # turn off legend name for ground_truth
  #     values = c(
  #       "TRUE" = "red",
  #       "FALSE" = "transparent")
  #   ) +
  #   scale_fill_gradient(low ="white",high =  "darkgreen")
  # 
  # make_escheR(sub_spe) |> 
  #   add_fill(var = "sum_umi_z", point_size=1.25) +
  #   # add_ground(var = "local_outliers",point_size=1.25, stroke = 1) +
  #   # scale_color_manual(
  #   #   name = "", # turn off legend name for ground_truth
  #   #   values = c(
  #   #     "TRUE" = "red",
  #   #     "FALSE" = "transparent")
  #   # ) +
  #   scale_fill_gradient(low ="white",high =  "darkgreen")
  
  
  sub_spe <- findArtifacts(
    sub_spe,
    mito_percent="expr_chrM_ratio",
    mito_sum="expr_chrM",
    n_rings=5,
    name="artifact"
  )
  
  print(
    make_escheR(sub_spe) |> 
    add_fill(var = "expr_chrM_ratio", point_size = 1.25) |> 
    add_ground(var = "artifact", point_size=1.25, stroke = 0.3) +
    scale_fill_gradient(low = "white", high = "black") +
    scale_color_manual(
      name = "Affected area", # turn off legend name for ground_truth
      values = c(
        "TRUE" = "red",
        "FALSE" = "transparent")
    ) +
    labs(title = .sample)
  
  )
  # Warning message:
  #   In .check_reddim_names(x, value, withDimnames) :
  #   non-NULL 'rownames(value)' should be the same as 'colnames(x)' for
  # 'reducedDim<-'. This will be an error in the next release of Bioconductor.
  
  
}
dev.off()



library(rbenchmark)
benchmark(
  `colData` = {colData(spe)[["sample_id"]]},
  `dollar` = spe$sample_id,
  replications = 100
)



# Session Info ------------------------------------------------------------
sessionInfo()

