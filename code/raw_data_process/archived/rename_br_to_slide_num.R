library(here)
library(tidyverse)

meta_df <- read_csv(
  file = here("code", "raw_data_process",
              "sample_meta_path.csv"))

expr_meta|> 
  pwalk(.f = function(BrNumbr, `Slide #`, `Array #`,
                      sr_fldr_path, ...){
    
    br_num_path <- file.path(processed_sparang_fldr, 
                             paste0("Br", BrNumbr, "_", `Array #`))
    
    if(!file.exists(br_num_path)){
      warning("No folder:", br_num_path)
      return()
    }
    
    rename_command <-
      paste("mv",
            br_num_path,
            sr_fldr_path
      )
    
    system(rename_command)

  })
