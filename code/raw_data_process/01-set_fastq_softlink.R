# Helper Function ---------------------------------------------------------
check_fastq_files_exist <- function(
    fastq_name_start,    # The start of fastq file name
    fastq_fldr_path      # The sample-specific fastq folder
) {
  tmp_fastq_files <- grep( 
    pattern = paste0(fastq_name_start,
                     "[^/]*",
                     "\\.fastq\\.gz$"),
    x = list.files(fastq_fldr_path),
    value = TRUE,
    ignore.case = TRUE # 2 samples have SHK instead of shk
  )
  
  return(length(tmp_fastq_files)==2)
}


# Create Sample Specific Folder -------------------------------------------
# Create "FASTQ" folder if not exists
mkdir_if_not_exist(
  processed_fastq_fldr, 
  recursive = TRUE # Create "raw-data" folder if not exists
)

# Error Prevention
stopifnot(dir.exists(processed_fastq_fldr))



# Create A Softlink -------------------------------------------------------
raw_all_fastq_files <- list.files(input_data_wh_folder, full.names = TRUE)

expr_meta |> 
  pwalk(.f = function(fastq_fldr_path,
                      fastq_name_start, ...){
    
    # Create Sample Directory
    mkdir_if_not_exist(fastq_fldr_path)
    
    # Check if the soft link is already set up
    if(check_fastq_files_exist(fastq_name_start, fastq_fldr_path))
      return() # Do nothing
    
    raw_fast_paths <- grep(
      pattern = paste0("/", fastq_name_start,
                       "[^/]*",
                       "\\.fastq\\.gz$"),
      x = raw_all_fastq_files,
      value = TRUE,
    )
    
    # Error Prevention
    # TODO: verify if it is true, there are 2 files
    stopifnot(length(raw_fast_paths)==2)
    
    # Create the command
    soft_copy_command <-
      paste("ln", "-s",
            paste(raw_fast_paths, collapse = " "),
            fastq_fldr_path
      )
    
    # Run the command
    system(soft_copy_command)
    
    # Validation step    
    stopifnot(
      check_fastq_files_exist(fastq_name_start, fastq_fldr_path)
    )
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



