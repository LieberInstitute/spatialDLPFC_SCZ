setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
  library(dplyr)
  library(tidyr)
})


############# neuropil ######################
#neuropil_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "neuropil_pseudo.rds"))
#
#dx_res <- registration_stats_enrichment(
#  neuropil_pseudo,
#  block_cor = NaN,
#  covars = c("age", "sex"),
#  var_registration = "neuropil_pos",
#  gene_ensembl = "gene_id",
#  gene_name = "gene_name"
#)
#
#write.csv(dx_res, file = here("processed-data", "image_processing", "enrichment", "neuropil_dx_res.csv"), row.names = TRUE)
#dx_res_sorted <- dx_res %>%
#  arrange(fdr_TRUE) %>%
#  filter(fdr_TRUE < 0.05 & logFC_TRUE>0)
#
#  dx_res_sorted <- dx_res %>%
#    arrange(fdr_TRUE) %>%
#    filter(p_value_TRUE < 0.05)
#head(dx_res_sorted)
spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
#
#### andrew code ####
###subset to 
brains = c("Br2719", "Br5182", "Br5367", "Br5472")
spe_sub = spe_ntc[, spe_ntc$brnum %in% brains]
mod = model.matrix(~neuropil_pos + PRECAST_07 + age + slide_id ,data =colData(spe_sub))

logcounts_matrix <- as.matrix(assays(spe_sub)$logcounts)
aligned_cells <- intersect(colnames(logcounts_matrix), rownames(mod))
logcounts_matrix <- logcounts_matrix[, aligned_cells, drop = FALSE]
mod <- mod[aligned_cells, , drop = FALSE]
fit <- bumphunter:::.getEstimate(logcounts_matrix, mod, coef=2, full=TRUE)

#fit = bumphunter:::.getEstimate(assays(spe_sub)$logcounts, mod, coef=2, full=TRUE)
tt = bumphunter:::.getModT(fit)
pv = 2*pt(-abs(tt$t), df = tt$df.total)
stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
rownames(stats) = rownames(assays(spe_sub)$logcounts)
stats$Symbol = rowData(spe_sub)$gene_name
save(stats, file = here("processed-data", "image_processing", "enrichment", "neuropil_dx_res_andrew_sub.Rdata"))

############# neuron ######################

#neun_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "neun_pseudo.rds"))
#dx_res <- registration_stats_enrichment(
#  neun_pseudo,
#  block_cor =NaN,
#  covars = c("age", "sex"),
#  var_registration = "neun_pos",
#  gene_ensembl = "gene_id",
#  gene_name = "gene_name"
#)
#
#write.csv(dx_res, file = here("processed-data", "image_processing", "enrichment", "neun_dx_res.csv"), row.names = TRUE)
#dx_res_sorted <- dx_res %>%
#  arrange(fdr_TRUE) %>%
#  filter(p_value_TRUE < 0.05)
#
#head(dx_res_sorted)
#
#spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
#colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
#
#### andrew code ####
###subset to 
#brains = c("Br2719", "Br5182", "Br5367", "Br5472")
#spe_sub = spe_ntc[, spe_ntc$brnum %in% brains]
#mod = model.matrix(~neun_pos + PRECAST_07 + age + slide_id ,data =colData(spe_sub))
#
#logcounts_matrix <- as.matrix(assays(spe_sub)$logcounts)
#aligned_cells <- intersect(colnames(logcounts_matrix), rownames(mod))
#logcounts_matrix <- logcounts_matrix[, aligned_cells, drop = FALSE]
#mod <- mod[aligned_cells, , drop = FALSE]
#fit <- bumphunter:::.getEstimate(logcounts_matrix, mod, coef=2, full=TRUE)
#
##fit = bumphunter:::.getEstimate(assays(spe_sub)$logcounts, mod, coef=2, full=TRUE)
#tt = bumphunter:::.getModT(fit)
#pv = 2*pt(-abs(tt$t), df = tt$df.total)
#stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
#rownames(stats) = rownames(assays(spe_sub)$logcounts)
#stats$Symbol = rowData(spe_sub)$gene_name
#save(stats, file = here("processed-data", "image_processing", "enrichment", "neun_dx_res_andrew_sub.Rdata"))
#
#
############## pnn ######################
#
#pnn_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "pnn_pseudo.rds"))
#dx_res <- registration_stats_enrichment(
#  pnn_pseudo,
#  block_cor =NaN,
#  covars = c("age", "sex"),
#  var_registration = "pnn_pos",
#  gene_ensembl = "gene_id",
#  gene_name = "gene_name"
#)
#
#write.csv(dx_res, file = here("processed-data", "image_processing", "enrichment", "pnn_dx_res.csv"), row.names = TRUE)
#dx_res_sorted <- dx_res %>%
#  arrange(fdr_TRUE) %>%
#  filter(p_value_TRUE < 0.05)
#
#head(dx_res_sorted)
#
##spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
##colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
##
##### andrew code ####
####subset to 
##brains = c("Br2719", "Br5182", "Br5367", "Br5472")
##spe_sub = spe_ntc[, spe_ntc$brnum %in% brains]
#mod = model.matrix(~pnn_pos + PRECAST_07 + age + slide_id ,data =colData(spe_sub))
#
#logcounts_matrix <- as.matrix(assays(spe_sub)$logcounts)
#aligned_cells <- intersect(colnames(logcounts_matrix), rownames(mod))
#logcounts_matrix <- logcounts_matrix[, aligned_cells, drop = FALSE]
#mod <- mod[aligned_cells, , drop = FALSE]
#fit <- bumphunter:::.getEstimate(logcounts_matrix, mod, coef=2, full=TRUE)
#
##fit = bumphunter:::.getEstimate(assays(spe_sub)$logcounts, mod, coef=2, full=TRUE)
#tt = bumphunter:::.getModT(fit)
#pv = 2*pt(-abs(tt$t), df = tt$df.total)
#stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
#rownames(stats) = rownames(assays(spe_sub)$logcounts)
#stats$Symbol = rowData(spe_sub)$gene_name
#save(stats, file = here("processed-data", "image_processing", "enrichment", "pnn_dx_res_andrew_sub.Rdata"))
#
############## claudin ######################
#
#vasc_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "vasc_pseudo.rds"))
#dx_res <- registration_stats_enrichment(
#  vasc_pseudo,
#  block_cor =NaN,
#  covars = c("age", "sex"),
#  var_registration = "vasc_pos",
#  gene_ensembl = "gene_id",
#  gene_name = "gene_name"
#)
#
#write.csv(dx_res, file = here("processed-data", "image_processing", "enrichment", "vasc_dx_res.csv"), row.names = TRUE)
#dx_res_sorted <- dx_res %>%
#  arrange(fdr_TRUE) %>%
#  filter(p_value_TRUE < 0.05)
#
#head(dx_res_sorted)
#
##spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
##colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)
##
##### andrew code ####
####subset to 
##brains = c("Br2719", "Br5182", "Br5367", "Br5472")
##spe_sub = spe_ntc[, spe_ntc$brnum %in% brains]
#mod = model.matrix(~vasc_pos + PRECAST_07 + age + slide_id ,data =colData(spe_sub))
#
#logcounts_matrix <- as.matrix(assays(spe_sub)$logcounts)
#aligned_cells <- intersect(colnames(logcounts_matrix), rownames(mod))
#logcounts_matrix <- logcounts_matrix[, aligned_cells, drop = FALSE]
#mod <- mod[aligned_cells, , drop = FALSE]
#fit <- bumphunter:::.getEstimate(logcounts_matrix, mod, coef=2, full=TRUE)
#
##fit = bumphunter:::.getEstimate(assays(spe_sub)$logcounts, mod, coef=2, full=TRUE)
#tt = bumphunter:::.getModT(fit)
#pv = 2*pt(-abs(tt$t), df = tt$df.total)
#stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
#rownames(stats) = rownames(assays(spe_sub)$logcounts)
#stats$Symbol = rowData(spe_sub)$gene_name
#save(stats, file = here("processed-data", "image_processing", "enrichment", "vasc_dx_res_andrew_sub.Rdata"))
