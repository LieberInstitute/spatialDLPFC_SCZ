### Prep scripts to check if VistoSeg::countnuclei returns the same counts from
###.   CN (the ones with centroid counts) mode and N (pixel based count)

library(here)
library(tidyverse)

samples <- list.files(here("processed-data/spaceranger")) |> 
  grep("^V[0-9]+F[0-9]+", x = _, value = TRUE)


count_path <- here("processed-data/spaceranger", samples,
                   "outs/spatial/tissue_spot_counts.csv") |> 
  set_names(samples)



count_path |> 
  map(
    .f = function(path){
      
      
      channel_names <- c("AF", "Claudin5", "DAPI", "NeuN", "WFA")
      
      
      if(!file.exists(path)){
        return(data.frame("AF" = NA, "Claudin5" = NA, "DAPI" = NA, "NeuN" = NA, "WFA" = NA))
      }
      
      count_df <- read.csv(path)
      
      channel_names |> set_names(channel_names) |> 
        map(.f = function(chnl){
          identical(
            getElement(count_df, paste0("N", chnl)),
            getElement(count_df, paste0("CN", chnl))
          )
        }) |> c() |> data.frame()
      
    }
  ) |>
  list_rbind(names_to = "Sample")

