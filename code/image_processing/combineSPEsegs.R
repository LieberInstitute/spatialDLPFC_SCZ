
spe = read.csv('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-053_A1/outs/spatial/tissue_spot_counts.csv')
spds = readRDS('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/rds/spatial_cluster/PRECAST/test_clus_label_df_semi_inform_k_2-16.rds')

library(tidyr)
spds$sample <- sapply(strsplit(spds$key, "-1_"), function(x) x[2])
spd = spds[spds$sample == "V12F14-053_A1", ]
spe$key = paste0(spe$barcode,"_","V12F14-053_A1")

sp <- merge(spe, spd, by = "key")


library("escheR")
library("ggplot2")
library("gridExtra")



p <- make_escheR(sp)
p1 <- p |> add_ground(var = "PRECAST_10")+scale_colour_viridis_d(option = "H")
p1 |> add_fill(var = "IWFA")+scale_fill_manual("greys")
print(p1)


ggplot(data = sp, aes(y = row, x = col, color = PRECAST_10, fill = IWFA))+geom_point(shape = )