library(here)
library(tidyverse)
library(cellranger)
library(glue)

# Load Helper Function ----------------------------------------------------
source(
  grep(
    pattern = "/fun[^/]*\\.R",
    x = list.files(
      here("code", "raw_data_process"),
      full.names = TRUE
    ),
    value = TRUE
  )
)

# Paths -------------------------------------------------------------------

exp_init <- "shk" |> tolower() # Seems always have lower case

# TODO: Scale to different round of experiment.
# TODO: Update this part
master_file_path <- here("raw-data", "sample_info",
                         "Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx") 

raw_file_path <- paste0(
  "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/",
  "raw-data/", 
  c("2023-01-29_SPag010323") #TODO: add here when more folder
)


fastq_fldr_path <- here("raw-data", "FASTQ")

processed_fldr_path <- here("processed-data", "spaceranger" )

# Import Sample Meta Data From Master Excel File --------------------------

meta_df <- readxl::read_excel(
  path = master_file_path,
  # "~/Downloads/Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx",
  col_names = TRUE,
  range = cellranger::cell_rows(
    c(2, NA) # TODO: check range
  ) # Start from the second row
)

# Simple error prevention
# Note: all `Sample #` should have value
stopifnot(all(!is.null(meta_df$`Sample #`)))


# Create file system data frame for folder names + fastq file names
fs_df <- meta_df |> 
  transmute(
    sample_fld_name = glue("Br{BrNumbr}_{`Array #`}"),
    raw_data_path = paste0(fastq_fldr_path, "/", sample_fld_name),
    process_data_path = paste0(processed_fldr_path, "/", sample_fld_name),
    # TODO: update "v" to match with 
    fastq_name_start = glue("{`Sample #`}v-{exp_init}"),
    # Information Needs for spaceranger script
    `Slide SN #`, `Array #`
  )


# Create Sample Specific Folder -------------------------------------------

# Create "FASTQ" folder if not exists

mkdir_if_not_exist(
  fastq_fldr_path, 
  recursive = TRUE # Create "raw-data" folder if not exists
)

# Error Prevention
stopifnot(dir.exists(fastq_fldr_path))

mkdir_if_not_exist(
  processed_fldr_path, 
  recursive = TRUE # Create "raw-data" folder if not exists
)

# Error Prevention
stopifnot(dir.exists(processed_fldr_path))


# Create A Softlink -------------------------------------------------------
# TODO: if there is multiple samples
all_files <- list.files(raw_file_path, full.names = TRUE)

fs_df |> 
  # pull(fastq_name_start) |> 
  pwalk(.f = function(raw_data_path,
                      fastq_name_start, ...){
    
    # Create Sample Directory
    mkdir_if_not_exist(raw_data_path)
    
    # Check if there the soft link set up
    tmp_files <- grep( 
      pattern = paste0("^", fastq_name_start,
                       "[^/]*",
                       "\\.fastq\\.gz$"),
      x = list.files(raw_data_path),
      value = TRUE,
      ignore.case = TRUE # 2 samples have SHK instead of shk
    )
    
    # .fastq files exists
    if(length(tmp_files)==2) return() # Do nothing
    
    raw_fast_paths <- grep(
      pattern = paste0("/", fastq_name_start,
                       "[^/]*",
                       "\\.fastq\\.gz$"),
      x = all_files,
      value = TRUE,
    )
    
    # Error Prevention
    # TODO: verify if it is true, there are 2 files
    stopifnot(length(raw_fast_paths)==2)
    
    # Create the command
    soft_copy_command <-
      paste("ln", "-s",
            paste(raw_fast_paths, collapse = " "),
            raw_data_path
      )
    
    # Run the command
    system(soft_copy_command)
    
    # Validation step    
    tmp_files <-grep( 
      pattern = paste0("^", fastq_name_start,
                       "[^/]*",
                       "\\.fastq\\.gz$"),
      x = list.files(raw_data_path),
      value = TRUE)
    
    stopifnot(length(tmp_files)==2)
  })



# Update folder permission
system(
  paste("sh",
        "/dcs04/lieber/lcolladotor/_jhpce_org_LIBD001/update_permissions_spatialteam.sh",
        fastq_fldr_path,
        sep = " ")
)



