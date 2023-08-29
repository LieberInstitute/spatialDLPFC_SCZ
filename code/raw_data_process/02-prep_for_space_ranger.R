# TODO: move this to space ranger part
# # Create Space Ranger folder if not exists
# mkdir_if_not_exist(
#   processed_sparang_fldr, 
#   recursive = TRUE # Create "raw-data" folder if not exists
# )

# stopifnot(dir.exists(processed_sparang_fldr))


# expr_meta

# Check success run -------------------------------------------------------

have_SR_summary <- function(smp_slide_id){
  file.exists(
    here(
      processed_sparang_fldr,
      smp_slide_id,
      "outs","web_summary.html")
  )
}

SR_folders_df <-
  data.frame(
    sample_name = grep(
      pattern = "^V",
      x = list.files(
        processed_sparang_fldr
      ),
      value = TRUE
    )
  ) |> 
  mutate(
    success = have_SR_summary(sample_name)
  ) |> 
  filter(success == TRUE)




# loupe_path <- here("processed-data", "VistoSeg", "loupe")

to_run_df <- dplyr::anti_join(
  x = expr_meta,
  y = SR_folders_df,
  by = "sample_name"
)

# save_path <- here("code", "spaceranger",
#  "parameters")
# mkdir_if_not_exist(save_path, recursive = TRUE)

# (Optional) TODO: clean up to run log and parameter files

#  Set up jobs ------------------------------------------------------------
to_run_df |> 
  pwalk(
    .f = function(
    sample_name,
    `Slide #`,
    `Array #`,
    loupe_file_path,
    fastq_fldr_path,
    ...
    ){
      browser()
      # for(i in 1:nrow(to_run_df)){
      # sample_name <- to_run_df$sample_name[i]
      
      # to_run_df |> 
      #   filter(sample_name == sample_name) |> 
      c(
        sample_name,
        `Slide #`,
        `Array #`,
        paste0(loupe_file_path, ".tif"),
        paste0(loupe_file_path, ".json"),
        fastq_fldr_path
      ) |> 
        write.table(
          file=file.path(pmtr_sparanger,
                         paste0(sample_name,".tsv")),
          quote=F, col.names=F, row.names=F, sep="\t",
          append = FALSE
        )
      
      # Assemble the job submission command
      
      # SGE version (Deprecated)
      # job_sub_commond <- paste(
      #   "qsub",
      #   "-N", paste0("run_SR_", sample_name), #Sample specific job name
      #   "-wd", lib_sparanger, # Starting directory
      #   file.path(lib_sparanger, "spaceranger_SGE.sh"),
      #   sep = " "
      # )
      
      
      # SLURM Implementation
      # TODO: write a SLURM version
      stop("SLURM version has not been implemented yet.")
      job_sub_commond <- paste(
        "qsub",
        "-N", paste0("run_SR_", sample_name), #Sample specific job name
        "-wd", lib_sparanger, # Starting directory
        file.path(lib_sparanger, "spaceranger_SLURM.sh"),
        sep = " "
      )
      
      
      # Run Spaceranger for each sample as a job
      system(
        job_sub_commond
      )
    }
    # }
  )




