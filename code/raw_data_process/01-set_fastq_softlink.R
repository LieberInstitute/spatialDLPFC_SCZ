# Helper Function (deprecated) ---------------------------------------------------------
# check_fastq_files_exist <- function(
#     fastq_name_start,    # The start of fastq file name
#     fastq_fldr_path      # The sample-specific fastq folder
# ) {
#   tmp_fastq_files <- grep( 
#     pattern = paste0(fastq_name_start,
#                      "[^/]*",
#                      "\\.fastq\\.gz$"),
#     x = list.files(fastq_fldr_path),
#     value = TRUE,
#     ignore.case = TRUE # 2 samples have SHK instead of shk
#   )
#   
#   return(length(tmp_fastq_files)==2)
# }


# Create Sample Specific Folder -------------------------------------------
# Create "FASTQ" folder if not exists
mkdir_if_not_exist(
  processed_fastq_fldr, 
  recursive = TRUE # Create "raw-data" folder if not exists
)

# Error Prevention
stopifnot(dir.exists(processed_fastq_fldr))



# Create A Softlink -------------------------------------------------------
# TODO: add new path to input_data_wh_folder in "file_pahts.R" when
#       new data comes from seq core
raw_all_fastq_files <- list.files(
  input_data_wh_folder,
  full.names = TRUE, recursive = TRUE)


expr_meta |> 
  pwalk(.f = function(
    fastq_fldr_path,
    # fastq_name_start,
    Experimenter,
    sample_name,
    `Sample #`,
    ...){
    # Create Sample Directory
    mkdir_if_not_exist(fastq_fldr_path)
    
    # Skip samples with softlink set up
    if(length(list.files(fastq_fldr_path)) != 0){
      return()
    }
    
    raw_fast_paths <- grep(
      # TODO: ask someone to clean up this regex
      pattern = paste0(
        "[^/]*", # File path prefix
        "/", `Sample #`, # Starts with /xxv or /PNN-xxv
        # "(/|`/PNN-`)", `Sample #`, # Starts with /xxv or /PNN-xxv
        "[-]", Experimenter, #e.g. `-shk` or `-SHK`
        "[^/]*",
        "\\.fastq\\.gz$"),
      x = raw_all_fastq_files,
      ignore.case = TRUE,
      value = TRUE
    )
    
    
    # Error Prevention
    
    # Exceptions: shk for dACC project
    # SHK did Visium experiment for dACC project too.
    if(any(str_detect(raw_fast_paths, "dACC"))){
      stop("Detected dACC sample in FASTQ files")
    }
    
    if(length(raw_fast_paths)==0){
      stop("No matched FASTQ file for Sample ", sample_name, 
           " with Experiment # ", `Sample #`)
    }
    
    # Create the command
    soft_copy_command <-
      paste("ln", "-s",
            paste(raw_fast_paths, collapse = " "),
            fastq_fldr_path
      )
    
    # Run the command
    system(soft_copy_command)
    
    # Validation step    
    # if there's at least 1 FASTQ in the folder
    if(length(list.files(fastq_fldr_path)) == 0){
      stop("Failed to set up soft link for Sample ", sample_name)
    }
  })


# Update folder permission
system(
  paste("sh",
        file.path(
          "/dcs04/lieber/lcolladotor/_jhpce_org_LIBD001",
          "update_permissions_spatialteam.sh"),
        processed_fastq_fldr,
        sep = " ")
)





