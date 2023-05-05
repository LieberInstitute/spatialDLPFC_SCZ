# Input File Path -------------------------------------

## Experiment Master File --------------------------------------------------
input_exp_raw_file <- here("raw-data", "experiment_info",
                         "VisiumSPG_PNN_Master.xlsx") 




## Data Warehouse Folder Path -----------------------------------------------
input_data_wh_folder <- file.path(
  # TODO: change drive path (if necessary) 
  "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001",
  "raw-data",  
  c(
    #TODO: add here when more subfolder
    "2023-01-29_SPag010323",
    "2022-12-13_SPag120622_2"
  ) 
)



# Outputs File Path -------------------------------------
# Default Path for Preprocess Outputs -------------------------------------
# No need to change

## fastq soft links  
processed_fastq_fldr <- here("raw-data", "FASTQ")

## spaceranger folder 
# TODO: update the name
processed_sparang_fldr <- here("processed-data", "spaceranger")
