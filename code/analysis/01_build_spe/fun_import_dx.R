# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(cellranger) # Extract excel file

source(here("code", "visium_data_process", "file_paths.R"))


# Import Excel File -------------------------------------------------------
raw_df <- readxl::read_excel(
  path = input_raw_dx_file,
  sheet = "Big_240_DLPFC_Dissections",
  col_names = TRUE,
  # First row is grouped cell names, hence starting second row
  range = cellranger::cell_rows(
    c(2, NA)
  ) # Start from the second row
) |> 
  # Only include relevent columns
  select(
    brain_num = `Brain #`,
    SCZ_study = `Study 1`,
    dx_raw = "Diagnosis", 
    age = "AGE",
    sex = "SEX",
    "PMI", "RIN",
    "Dissection Date"
    #  Seems to be notes
    # "opiate dx"
  ) 
# Remove example row based on "Ex" --------------------------------------
ex_row_id <- grep(
  pattern = "^Ex:.*",
  x = raw_df$brain_num)

clean_df <- raw_df[-ex_row_id,] |> 
  filter_all(any_vars(!is.na(.))) # Remove empty rows.

# Validation
# if row starts with "Ex:" in dx column
stopifnot(
  grep(
    pattern = "^Ex:.*",
    x = clean_df$dx_raw) |> 
    length() == 0
)

# Clean up brain number column--------------------------------------------------
# Remove all " - xxNotesxx"
clean_df$brain_num  <- str_sub(
  clean_df$brain_num,
  1, 6
)


# Validation
# If all brain numbers are 6 characters
# This validation is not robust to longer Brain nums
stopifnot(
  all(str_length(clean_df$brain_num)==6)
)

# Validation

clean_df <- clean_df[!duplicated(clean_df$brain_num), ]

stopifnot(all(!duplicated(clean_df$brain_num)))
# clean_df |> group_by(brain_num) |> summarize(n = n()) |> filter(n!=1)

# Clean up dx column --------------------------------------------------

# Helper functions
# Change all elements starting with "sc" to scz (lower case)
clean_scz <- function(vec){
  scz_ind <- grep(
    pattern = "^sc.+",
    x = vec
  )
  
  vec[scz_ind] <- "scz"
  
  return(vec)
}

clean_ntc <- function(vec){
  browser()
  ntc_ind <- str_detect(
    string =  vec,
    pattern = "control|neurotypical|ntc"
  )
  
  vec[ntc_ind] <- "ntc"
  
  return(vec)
}

clean_df <- clean_df |> 
  mutate(
    dx = case_when(
      str_detect(tolower(dx_raw), "^sc.+") ~ "scz",
      str_detect(tolower(dx_raw), "control|neurotypical|ntc") ~ "ntc",
      # Set other dx to NA
      TRUE ~ NA_character_
    ),
    race = "CAUC"  # NOTE: "all donors in this study are European ancestry"
  ) |> 
  filter(!is.na(dx)) |> 
  select(-c(SCZ_study, dx_raw))



