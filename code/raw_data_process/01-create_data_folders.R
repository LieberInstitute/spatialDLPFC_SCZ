library(here)
library(tidyverse)
library(cellranger)
library(glue)


# Import Sample Meta Data From Master Excel File --------------------------
exp_init <- "shk"


# TODO: Scale to different round of experiment.
master_file_path <- here("raw-data", "sample_info",
                         "Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx") 

raw_file_path <- paste0("/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/",
                        "raw-data/", "2023-01-29_SPag010323")

meta_df <- readxl::read_excel(
  path = master_file_path,
  # "~/Downloads/Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx",
  col_names = TRUE,
  range = cellranger::cell_rows(
    c(2, NA)
  ) # Start from the second row
)

# Simple error prevention
# Note: all `Sample #` should have value
stopifnot(all(!is.null(meta_df$`Sample #`)))


# Create file system data frame for folder names + fastq file names
fs_df <- meta_df |> 
  transmute(
    sample_fld_name = glue("Br{BrNumbr}_{`Array #`}"),
    # TODO: update "v" to match with 
    fastq_name_start = glue("{`Sample #`}v-{exp_init}")
  )


# Create Sample Specific Folder -------------------------------------------

fastq_fldr_path <- here("raw-data", "FASTQ")

# Create "FASTQ" folder if not exists

if(!dir.exists(fastq_fldr_path))
  dir.create(
    fastq_fldr_path, 
    recursive = TRUE # Create "raw-data" folder if not exists
  )

# Error Prevention
stopifnot(dir.exists(fastq_fldr_path))


# Create Sample Specific Folder
fs_df |> 
  pull(sample_fld_name) |> 
  walk(
    ~mkdir_if_not_exist(
      here(fastq_fldr_path, .x)
    )
  )


# Create A Softlink -------------------------------------------------------
# TODO: if there is multiple samples
all_files <- list.files(raw_file_path, full.names = TRUE)

fs_df |> 
  pull(fastq_name_start) |> 
  walk(.f = function(name_start){
    raw_fast_paths <- grep(
      pattern = paste0("^", name_start,
                       ".*",
                       "\.fastq.gz"),
      x = name_start,
      value = TURE,
    )
    
    # TODO: verify if it is true, there are 2 files
    stopifnot(length(raw_fast_paths)!=2)
    
    # Create the command
    soft_copy_command <-
      paste("ln", "-s",
            paste(raw_fast_path, collapse = " "),
            here(fastq_fldr_path, .x)
            )
    
    # Run the command
    system(soft_copy_command)
    
    # TODO:Validation step
    
  })




# Update folder permission
system(
  paste("sh",
        "/dcs04/lieber/lcolladotor/_jhpce_org_LIBD001/update_permissions_spatialteam.sh",
        fastq_fldr_path,
        sep = " ")
)



