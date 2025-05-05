setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
scz = readRDS(here("processed-data", "image_processing", "enrichment", "spe_scz.rds"))

df_ntc = as.data.frame(colData(ntc))
df_scz = as.data.frame(colData(scz))

dfN = df_ntc %>%
summarize(Neuropil_spots = sum(neuropil_pos == TRUE), NeuN_spots = sum(neun_pos == TRUE), 
WFA_spots = sum(pnn_pos == TRUE), Claudin5_spots = sum(vasc_pos == TRUE))

#dfN
#  Neuropil_spots NeuN_spots WFA_spots Claudin5_spots
#          71804      28957      7084          10708

dfS = df_scz %>%
summarize(Neuropil_spots = sum(neuropil_pos == TRUE), NeuN_spots = sum(neun_pos == TRUE), 
WFA_spots = sum(pnn_pos == TRUE), Claudin5_spots = sum(vasc_pos == TRUE))

#dfS
#  Neuropil_spots NeuN_spots WFA_spots Claudin5_spots
#          76161      26060      5249          11198

dfNT <- df_ntc %>%
group_by(PRECAST_07) %>%
summarize(NNeuropil = sum(neuropil_pos == TRUE), PNeuropil = NNeuropil/dfN$Neuropil_spots, NNeuN = sum(neun_pos == TRUE), PNeuN = NNeuN/dfN$NeuN_spots,
NWFA = sum(pnn_pos == TRUE), PWFA = NWFA/dfN$WFA_spots, NClaudin5 = sum(vasc_pos == TRUE), PClaudin5 = NClaudin5/dfN$Claudin5_spots)

dfSC <- df_scz %>%
group_by(PRECAST_07) %>%
summarize(NNeuropil = sum(neuropil_pos == TRUE), PNeuropil = NNeuropil/dfS$Neuropil_spots, NNeuN = sum(neun_pos == TRUE), PNeuN = NNeuN/dfS$NeuN_spots,
NWFA = sum(pnn_pos == TRUE), PWFA = NWFA/dfS$WFA_spots, NClaudin5 = sum(vasc_pos == TRUE), PClaudin5 = NClaudin5/dfS$Claudin5_spots)
