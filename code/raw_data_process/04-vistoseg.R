
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

sr_success_brs <- file.exists(
  here("processed-data", "spaceranger",
       finished_df$sample_fld_name,
       "outs","web_summary.html")
)

# TODO: it is possible to just run successful_runs
if(nrow(finished_df) != sum(sr_success_brs)) # Space ranger Fails exist
  stop(
    paste0("Failed space ranger runs for ",
           paste0(
             finished_df$sample_fld_name[!sr_success_brs],
             collapse = ", ")
    )
  )


# Locate DAPI Segmentation Files ------------------------------------------

# dapi_seg_data_files <- list.files(
#   path = here("processed-data", "RealPNN", "capture_area_segmentations", 
#               "DAPI", "Data_files")
# )

vistoseg_finished_brs <- file.exists(
  here("processed-data", "spaceranger",
       finished_df$sample_fld_name,
       "outs","spatial", "tissue_positions_list.csv")
)


# Find the Br that needs to run VistoSeg
# No vistoseg output + success sr run

dapi_fs_df <- finished_df[sr_success_brs & (!vistoseg_finished_brs), , drop = FALSE] |>  
  left_join(
    fs_df,
    by = "sample_fld_name"
  ) |> transmute(
    sample_name = sample_fld_name,
    sr_fld_path = process_data_path,
    dapi_file_name = glue("{`Slide SN #`}_{`Array #`}_info.csv"),
    dapi_file_path = here("processed-data", "RealPNN", "capture_area_segmentations", 
                          "DAPI", "Data_files", dapi_file_name)
  )


have_dapi_file <- dapi_fs_df |>
  pull(dapi_file_path)|> 
  file.exists()

if(nrow(dapi_fs_df) > sum(have_dapi_file)) {
  paste0(
    "Some samples doesn't have DAPI files, and hence not runing VistoSeg: ",
    paste0(dapi_fs_df$sample_name[!have_dapi_file],  collapse = ", ")
  ) |> 
    warning()
}


vistseg_fld_path <- "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg"

mkdir_if_not_exist(here(vistseg_fld_path, "parameters"))
mkdir_if_not_exist(here(vistseg_fld_path, "logs"))
                        

# Create Sample Specific VistoSeg Job -------------------------------------

dapi_fs_df[have_dapi_file, , drop = FALSE] |> 
  pwalk(.f = function(sample_name, dapi_file_path, sr_fld_path, ...){
    browser()
    # Create tsv file for sample parameters
    data.frame(
      neuc_seg_mat = dapi_file_path,
      json = here(sr_fld_path,
                  "outs/spatial/scalefactors_json.json"),
      output_path = here(sr_fld_path,
                         "outs/spatial/tissue_positions_list.csv")
    ) |> 
      write.table( 
        file=here(
          "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/vistoseg",
          "parameters",
          paste0(sample_name,".tsv")
        ),
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
  )
