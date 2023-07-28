library(SpatialExperiment)
library(tidyverse)
library(escheR)

spe <- readRDS(
  here("processed-data", "rds", "spe", "test_spe_bayesSpace_harmony.rds")
)


# Choose the SpD ----------------------------------------------------------

spe$MBP <- logcounts(spe)["ENSG00000197971",]
spe$spd <- factor(spe$bayesSpace_harmony_9 %in% c(2,3) , levels = c(TRUE))

unique(spe$sample_id) |> 
  map(
    .f = function(id){
      # browser()
      uni_spe <- spe[, spe$sample_id == id]
      stopifnot(ncol(uni_spe)!=0)
      
      # uni_spe$bayesSpace_harmony_9 <- factor(uni_spe$bayesSpace_harmony_9)
      # This is logcount
      # make_escheR(uni_spe) |> 
      #   add_fill(var = "bayesSpace_harmony_9")
      
      make_escheR(uni_spe) |> 
        add_fill(var = "MBP") |> 
        # TODO: check 
        add_ground(var = "spd") +
        scale_color_manual(values = c("red"),
                           na.value = "transparent")
    }
   
  )


