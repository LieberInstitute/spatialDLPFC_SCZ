setwd("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")
library(readxl)

Dr = here()
sang_ho = as.data.frame(read_excel(here("processed-data", "image_processing","DLPFC_240_collection_updated_final.xlsx")))
data = as.data.frame(read_excel(here("raw-data","experiment_info","VisiumSPG_PNN_Master.xlsx")))

data$Brnumbr = paste0("Br",data$BrNumbr)

key = match(sang_ho$Brnumbr, data$Brnumbr)

dataNew = data.frame(array = paste0(data[,6] , "_", data[,7]), brnum = sang_ho$Brnumbr, Dx = sang_ho$Assigned_Diagnosis)

write.csv(dataNew, here("processed-data", "image_processing", "Dxinfo.csv"))