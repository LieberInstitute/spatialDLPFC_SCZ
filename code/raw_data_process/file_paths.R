# Experiment Master File --------------------------------------------------

master_file_path <- here("raw-data", "sample_info",
                         "Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx") 




# Data Warehouse Folder Path -------------------------------------------------------
raw_folder_path <- file.path(
  # TODO: change drive path (if necessary) 
  "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001",
  "raw-data",  
  c(
    #TODO: add here when more subfolder
    "2023-01-29_SPag010323",
    "2022-12-13_SPag120622_2"
  ) 
)


# Default Path for Preprocess Outputs -------------------------------------
# No need to change

# fastq folder saving the soft links
fastq_fldr_path <- here("raw-data", "FASTQ")

# spaceranger folder
# TODO: update the name
processed_fldr_path <- here("processed-data", "spaceranger")
