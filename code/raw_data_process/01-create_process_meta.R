# Import Sample Meta Data From Master Excel File --------------------------
raw_df_expr_meta <- readxl::read_excel(
  path = input_exp_raw_file,
  col_names = TRUE,
  sheet = "Summary"
)

# Error prevention
# Note: all `Sample #` should have value
stopifnot(all(!is.null(raw_df_expr_meta$`Sample #`)))


expr_meta <- raw_df_expr_meta |> 
  select(`Experiment #`:`Array #`) |>
  # rename_all(tolower) |> 
  mutate(
    BrNumbr = paste0("Br", BrNumbr),
    sample_name = glue("{`Slide #`}_{`Array #`}"),
    # Fastq Softlinks Folder
    fastq_fldr_path = file.path(processed_fastq_fldr, sample_name),
    # Space Ranger Folder
    sr_fldr_path = file.path(processed_sparang_fldr, sample_name),
    # Loupe Browser Folder (Pipe into Space Ranger)
    loupe_file_path = file.path(processed_loupe_fldr, sample_name)
    # TODO: add other
  )

write.csv(expr_meta,
          file = here("code", "raw_data_process",
                      "sample_meta_path.csv"),
          append = FALSE, row.names = FALSE, col.names = FALSE)