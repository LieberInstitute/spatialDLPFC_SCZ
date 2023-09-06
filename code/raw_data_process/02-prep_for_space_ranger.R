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

# expr_meta <- readr::read_csv(
#   file = here("code", "raw_data_process",
#               "sample_meta_path.csv")
# )


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
  # filter(`Sample #` == "13v") |>  # TODO: delete for latter
  pwalk(
    .f = function(
    sample_name,
    `Slide #`,
    `Array #`,
    loupe_file_path,
    fastq_fldr_path,
    ...
    ){
      # browser()

      c(
        sample_name,
        `Slide #`,
        `Array #`,
        paste0(loupe_file_path, ".tif"),
        paste0(loupe_file_path, ".json"),
        fastq_fldr_path
      ) |>
        t() |>
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
      job_sub_commond <- paste(
        "sbatch",
        "--job-name", sample_name, #Sample specific job name
        "--chdir", lib_sparanger, # Starting directory
        file.path(lib_sparanger, "spaceranger_SLURM.sh"),
        sep = " "
      )
      
      stop("SLURM version has not been tested yet.")
      
      # Run Spaceranger for each sample as a job
      system(
        job_sub_commond
      )
    }
    # }
  )




