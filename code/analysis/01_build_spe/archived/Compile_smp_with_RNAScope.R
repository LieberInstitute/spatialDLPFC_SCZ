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
  filter(
    CALL %in% c("PASS", "PASS*")
  ) |> 
  # Only include relevent columns
  select(
    brain_num = `Brain #`,
    SCZ_study = `Study 1`,
    dx_raw = "Diagnosis", 
    age = "AGE",
    sex = "SEX",
    "PMI", "RIN",
    "Dissection Date",
    RNAScope
    #  Seems to be notes
    # "opiate dx"
  ) 
# Remove example row based on "Ex" --------------------------------------
# ex_row_id <- grep(
#   pattern = "^Ex:.*",
#   x = raw_df$brain_num)

# clean_df <- raw_df[-ex_row_id,] |> 
#   filter_all(any_vars(!is.na(.))) # Remove empty rows.

clean_df <- raw_df

# Validation
# if row starts with "Ex:" in dx column
stopifnot(
  grep(
    pattern = "^Ex:.*",
    x = clean_df$dx_raw) |> 
    length() == 0
)


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



# Clean up brain number column--------------------------------------------------
# Remove all " - xxNotesxx"
clean_df$brain_num_only  <- str_sub(
  clean_df$brain_num,
  1, 6
)


input_exp_raw_file <- here("raw-data", "experiment_info",
                           "VisiumSPG_PNN_Master.xlsx") 
study_df <- readxl::read_excel(
  path = input_exp_raw_file,
  col_names = TRUE,
  sheet = "Summary"
) |> 
  transmute(BrNumbr = paste0("Br", BrNumbr))


clean_df |> 
  filter(
    brain_num_only %in%study_df$BrNumbr,
    !is.na(RNAScope),
    brain_num_only!="Br3942",
    brain_num != "Br2719 - depleted") |> 
  arrange(brain_num) |> 
  group_by(dx) |> 
  summarize(n = n())


# dx        n
# <chr> <int>
#   1 ntc      12
#   2 scz       8


