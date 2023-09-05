source(here("code", "raw_data_process", "fun_mkdir_if_not_exist.R"))

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
    # NOTE: should not contain the dash at the end
    # Wrong example "2023-06-29_KMay061223/KMay061223-391095963/"
    # Correct example "2023-06-29_KMay061223/KMay061223-391095963"
    "2023-01-29_SPag010323",
    "2022-12-13_SPag120622_2",
    "2023-06-29_KMay061223/KMay061223-391095963",
    "2023-08-24_SPag080923_S4"
  ) 
)


# Dx File Path ---------------------------------------------------------------
input_raw_dx_file <- "~/DLPFC cross-disorders brain collection bookkeeping updated 03_28_2023.xlsx"


# Outputs File Path -------------------------------------
# Default Path for Preprocess Outputs -------------------------------------
# No need to change

## fastq soft links  
processed_fastq_fldr <- here("raw-data", "FASTQ")

## spaceranger folder 
processed_sparang_fldr <- here("processed-data", "spaceranger")

## Loupe Output folder 
processed_loupe_fldr <- here("processed-data", "VistoSeg", "loupe")



# Parameter Files ---------------------------------------------------------


## Lib Path ----------------------------------------------------------------
# Space Ranger Path
lib_sparanger <- here("code", "raw_data_process", "spaceranger")


pmtr_sparanger <- file.path(lib_sparanger, "parameters")
mkdir_if_not_exist(pmtr_sparanger)


#TODO: create these folders

