# TODO: move this to space ranger part
# # Create Space Ranger folder if not exists
# mkdir_if_not_exist(
#   processed_sparang_fldr, 
#   recursive = TRUE # Create "raw-data" folder if not exists
# )

# stopifnot(dir.exists(processed_sparang_fldr))


# See which samples needs to be run with 
finished_df <- 
  data.frame(
    sample_fld_name = grep(
      pattern = "^Br",
      x = list.files(
        here("processed-data", "spaceranger")
      ),
      value = TRUE
    )
  )



# Add error prevention to check if
# the spaceranger run is acutally finished

sr_success_brs <- file.exists(
  here("processed-data", "spaceranger",
       finished_df$sample_fld_name,
       "outs","web_summary.html")
)

finished_df <- finished_df[sr_success_brs, , drop = FALSE]

loupe_path <- here("processed-data", "VistoSeg", "loupe")

to_run_df <- dplyr::anti_join(
  x = fs_df,
  y = finished_df,
  by = "sample_fld_name"
) |> 
  transmute(
    sample_name  = sample_fld_name,
    slide_id =  `Slide SN #`,
    array_id = `Array #`,
    img_path = file.path(loupe_path, glue("{slide_id}_{array_id}.tif")),
    json_path =file.path(loupe_path, glue("{slide_id}_{array_id}.json")),
    fastq_path = raw_data_path,
  )

save_path <- here("code", "spaceranger",
                  "parameters")

mkdir_if_not_exist(save_path, recursive = TRUE)

if(nrow(to_run_df) != 0){
  for(i in 1:nrow(to_run_df)){
    sample_name <- to_run_df$sample_name[i]
    write.table(to_run_df[i,], 
                file=paste0(save_path,"/",
                            sample_name,".tsv"),
                quote=F, col.names=F, row.names=F, sep="\t",
    )
    
    # Assemble the job submission command
    job_sub_commond <- paste(
      "qsub",
      "-N", sample_name, #Sample specific job name
      "-wd", here("code", "spaceranger"), # Starting directory
      here("code", "spaceranger", "spaceranger.sh"),
      sep = " "
    )
    
    
    # Run Spaceranger for each sample as a job
    system(
      job_sub_commond
    )
  }
}



