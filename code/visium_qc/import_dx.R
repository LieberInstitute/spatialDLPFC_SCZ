# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(cellranger) # Extract excel file



# File Path ---------------------------------------------------------------
path_raw_dx <- "~/DLPFC cross-disorders brain collection bookkeeping updated 03_28_2023.xlsx"

# Import Excel File -------------------------------------------------------
raw_df <- readxl::read_excel(
  path = path_raw_dx,
  sheet = "Big_240_DLPFC_Dissections",
  col_names = TRUE,
  # First row is grouped cell names, hence starting second row
  range = cellranger::cell_rows(
    c(2, NA) # TODO: check range
  ) # Start from the second row
) |> 
  # Only include relevent columns
  select(
    brain_num = `Brain #`,
    dx_raw = "Diagnosis", 
    age = "AGE",
    sex = "SEX"#,
    # TODO: what variables these are
    # "PMI", "RIN",
    #  Seems to be notes
    # "opiate dx"
  )

# Remove example row based on "Ex" --------------------------------------
ex_row_id <- grep(
  pattern = "^Ex:.*",
  x = raw_df$brain_num)

clean_df <- raw_df[-ex_row_id,]

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
clean_df$brain_num  <-  gsub(
  pattern = " - .*$",
  replacement = "",
  x = clean_df$brain_num
)

# Validation
# If all brain numbers are 6 characters
# This validation is not robust to longer Brain nums
stopifnot(
  all(str_length(clean_df$brain_num)==6)
)


# Clean up dx column

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
    )
  ) |> 
  filter(!is.na(dx))

# table(clean_df$dx)


# Descriptive Table -------------------------------------------------------
library(gtsummary)

# expr_meta <- read_csv()


# TODO: merge with expr_meta brain informaiton
clean_df %>%
  select(age, sex, dx) |> 
  tbl_summary(by = dx)


