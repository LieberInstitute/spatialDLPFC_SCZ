library(STexampleData)
library(SpatialExperiment)
library(scuttle)
library(scater)

library(escheR)

spe <- Visium_humanDLPFC()
spe <- spe[, spe$in_tissue == 1]

spe <- logNormCounts(spe)

positiv

RBG_norm_MBP <- scales::rescale(
  # TODO: change this to logcounts
  logcounts(spe)[which(rowData(spe)$gene_name == "MBP"),] |>
   scale(center = median(_), scale = FALSE) |>
   sapply(max, 0),
   to = c(0, 1))
RBG_norm_PCP4 <- scales::rescale(
  logcounts(spe)[which(rowData(spe)$gene_name == "PCP4"),] |>
   scale(center = median(_), scale = FALSE) |>
   sapply(max, 0),
   to = c(0, 1))
RBG_norm_SNAP25 <- scales::rescale(
  logcounts(spe)[which(rowData(spe)$gene_name == "SNAP25"),] |>
   scale(center = median(_), scale = FALSE) |>
   sapply(max, 0),
   to = c(0, 1)
   )


spe$rbg_val <- rgb(
  red = RBG_norm_MBP, 
  green = RBG_norm_PCP4,
  blue = RBG_norm_SNAP25)


make_escheR(spe) |> 
add_fill("rbg_val") +
scale_fill_identity() +
theme(plot.background = element_rect(fill = "black"))


spe$logcounts_SNAP25 <- logcounts(spe)[which(rowData(spe)$gene_name == "SNAP25"),]

make_escheR(spe) |> 
add_fill("logcounts_SNAP25") 