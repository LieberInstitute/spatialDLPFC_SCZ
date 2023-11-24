library("here")
library("readxl")
library("sessioninfo")

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


if(nrow(finished_df) > sum(sr_success_brs)) # Fails exist
  stop(
    paste0("Failed space ranger runs for ",
            paste0(
              finished_df$sample_fld_name[!sr_success_brs],
              collapse = ", ")
            )
    )



# Space Ranger QC ---------------------------------------------------------
## Adapted from https://github.com/LieberInstitute/ranger_metrics/tree/master/code/spaceranger

master_file_path <- here("raw-data", "sample_info",
        "Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx")

# TODO: remove this part. This part is for test purpose
# info <- readxl::read_excel(
#   path = master_file_path,
#   # "~/Downloads/Visium_IF_SCZ_PNN_1stRound_08142022_MasterExcel.xlsx",
#   col_names = TRUE,
#   range = cellranger::cell_rows(
#     c(2, NA) # TODO: check range
#   ) # Start from the second row
# )

# TODO:

stop("not implemented yet")
# info <- readxl::read_excel(master_file_path, sheet = 1)
# sample_col <- grep("sample \\#", tolower(colnames(info)))
# stopifnot(length(sample_col) == 1)
# 
# n_sample <- sum(!is.na(info[, sample_col])) - sum(tolower(info[, sample_col]) == "sample")
# 
# info <-
#   info[, !grepl("dilution|^\\.+", tolower(colnames(info)))]
# 
