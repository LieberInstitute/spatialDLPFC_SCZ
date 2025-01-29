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

# Load Data ----
load(file = here("processed-data", "image_processing", "enrichment", "neuropil_dx_res_andrew_sub.Rdata"))
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "neuropil_dx_res.csv"))
neuropil = read_excel("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/gene_sets/aau3644_TableS1.xls", sheet="DE", skip=4)
neuropil = as.data.frame(neuropil)
neuropil$padj = as.numeric(neuropil$padj)

spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)

### Andrew stats and homer data
#table(neuropil$padj < 0.05, sign(neuropil$stat))
stats$S_logFC = neuropil$Log2FoldChange[match(stats$Symbol, toupper(neuropil$GeneID))]
stats$S_padj = neuropil$padj[match(stats$Symbol, toupper(neuropil$GeneID))]

ct = cor.test(stats$logFC, stats$S_logFC)
ct_sig = with(stats[which(stats$S_padj < 0.05),], cor.test(logFC, S_logFC))
ct_fil = with(stats[which(stats$P.Value < 0.05),], cor.test(logFC, S_logFC))
ct_both = with(stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),], cor.test(logFC, S_logFC))

n_all = sum(complete.cases(stats[,c("logFC", "S_logFC")]))
n_sig = sum(complete.cases(stats[which(stats$S_padj < 0.05),c("logFC", "S_logFC")]))
n_fil = sum(complete.cases(stats[which(stats$P.Value < 0.05),c("logFC", "S_logFC")]))
n_both = sum(complete.cases(stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),c("logFC", "S_logFC")]))

pdf(here("plots", "image_processing", "enrichment", "neuropil_andrewSub_homer.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)

plot(logFC ~ S_logFC, data=stats,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed homologs"),
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$S_padj < 0.05),],
	main = paste(n_sig, "Neuropil Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05),],
	main = paste(n_fil, "-DAPI Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),],
	main = paste(n_both, "-DAPI & Neuropil Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()


# ## by layer
# layerIdx = splitit(sce$layer_guess)
# statList_layer = mclapply(layerIdx,  function(ii) {
	# mod = model.matrix(~zero_cell + subject_position,
		# data =colData(sce[,ii]))
	# fit = bumphunter:::.getEstimate(
		# sce_logcounts[,ii], mod, coef=2, full=TRUE)
	# tt = bumphunter:::.getModT(fit)
	# pv = 2*pt(-abs(tt$t), df = tt$df.total)
	# stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
	# rownames(stats) = rownames(sce_logcounts)
	# stats
# }, mc.cores=7)

# save(stats, statList_layer,
	# file = "rda/linear_model_zerospot_byLayer.Rdata")
	

neuropil <- read_excel(here("raw-data/images/SPG_Spot_Valid/Neuropil/Supple_Table_10_Neuropil.xlsx"), sheet = "10_human_synase_markers", skip = 1)
# Remove the first 8 columns
neuropil <- neuropil[, -(1:8)]
neuropil = as.data.frame(neuropil)

## with andrew stats ##
stats$S_logFC = neuropil$avg_log2FC...11[match(stats$Symbol, neuropil$gene...9)]
stats$S_padj = neuropil$p_val_adj...14[match(stats$Symbol, neuropil$gene...9)]

ct = cor.test(stats$logFC, stats$S_logFC)
ct_sig = with(stats[which(stats$S_padj < 0.05),], cor.test(logFC, S_logFC))
ct_fil = with(stats[which(stats$P.Value < 0.05),], cor.test(logFC, S_logFC))
ct_both = with(stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),], cor.test(logFC, S_logFC))

n_all = sum(complete.cases(stats[,c("logFC", "S_logFC")]))
n_sig = sum(complete.cases(stats[which(stats$S_padj < 0.05),c("logFC", "S_logFC")]))
n_fil = sum(complete.cases(stats[which(stats$P.Value < 0.05),c("logFC", "S_logFC")]))
n_both = sum(complete.cases(stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),c("logFC", "S_logFC")]))

stats_sig = stats[which(stats$P.Value < 0.05),]
expr_data <- assay(spe_ntc, "logcounts")[rownames(stats_sig), ]
sum_expression <- colSums(expr_data)

pdf(here("plots", "image_processing", "enrichment", "neun_andrewSub_muchun.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = neuropil_pos, y = sum_expression)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA)
  labs(title = "Sum of Expression of -DAPI DEGs",
       x = "-DAPI",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(logFC ~ S_logFC, data=stats,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$S_padj < 0.05),],
	main = paste(n_sig, "Neuropil Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05),],
	main = paste(n_fil, "-DAPI Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),],
	main = paste(n_both, "Neuropil & -DAPI Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_stats = stats %>% filter(P.Value < 0.05) %>% arrange(P.Value)
head(filtered_stats,20)

          1 
#                     logFC         t       P.Value    Symbol   S_logFC
#ENSG00000198888  0.3770988  40.05265  0.000000e+00    MT-ND1        NA
#ENSG00000198763  0.3670460  40.04855  0.000000e+00    MT-ND2        NA
#ENSG00000198804  0.2781078  39.45632  0.000000e+00    MT-CO1        NA
#ENSG00000198712  0.3559335  40.15322  0.000000e+00    MT-CO2        NA
#ENSG00000198899  0.3491187  41.02401  0.000000e+00   MT-ATP6        NA
#ENSG00000198938  0.3261531  40.66989  0.000000e+00    MT-CO3        NA
#ENSG00000198840  0.3390771  42.63933  0.000000e+00    MT-ND3        NA
#ENSG00000198727  0.3410232  38.03770 4.813629e-303    MT-CYB        NA
#ENSG00000198886  0.3008395  37.50908 4.401499e-295    MT-ND4        NA
#ENSG00000198786  0.3571992  32.70693 4.107491e-227    MT-ND5        NA
#ENSG00000269028  0.3480850  31.06119 1.118704e-205 MTRNR2L12 1.4801416
#ENSG00000080824 -0.2761411 -24.18130 6.892873e-127  HSP90AA1 0.5682143
#ENSG00000228253  0.2934994  22.67738 4.520686e-112   MT-ATP8        NA
#ENSG00000212907  0.2669093  21.85816 2.252183e-104   MT-ND4L        NA
#ENSG00000251562 -0.3583321 -20.73987  2.696844e-94    MALAT1 1.2789470
#ENSG00000166598 -0.2440817 -19.00174  1.265148e-79   HSP90B1        NA
#ENSG00000127914 -0.2202216 -18.94875  3.390887e-79     AKAP9        NA
#ENSG00000167996  0.1682857  18.47478  2.039074e-75      FTH1 0.7537682
#ENSG00000089737 -0.2314099 -18.08225  2.345337e-72     DDX24 0.2723909
#ENSG00000085224 -0.2262659 -17.62009  7.824319e-69      ATRX 0.2708158

### Boyi stats 
dx_res$S_logFC = neuropil$avg_log2FC...11[match(dx_res$gene, neuropil$gene...9)]
dx_res$S_padj = neuropil$p_val_adj...14[match(dx_res$gene, neuropil$gene...9)]

ct = cor.test(dx_res$logFC_TRUE, dx_res$S_logFC)
ct_sig = with(dx_res[which(dx_res$S_padj < 0.05),], cor.test(logFC_TRUE, S_logFC))
ct_fil = with(dx_res[which(dx_res$p_value_TRUE < 0.05),], cor.test(logFC_TRUE, S_logFC))
ct_both = with(dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),], cor.test(logFC_TRUE, S_logFC))

n_all = sum(complete.cases(dx_res[,c("logFC_TRUE", "S_logFC")]))
n_sig = sum(complete.cases(dx_res[which(dx_res$S_padj < 0.05),c("logFC_TRUE", "S_logFC")]))
n_fil = sum(complete.cases(dx_res[which(dx_res$p_value_TRUE < 0.05),c("logFC_TRUE", "S_logFC")]))
n_both = sum(complete.cases(dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),c("logFC_TRUE", "S_logFC")]))

dx_res_sig = dx_res[which(dx_res$S_padj < 0.05),]
expr_data_sig <- assay(spe_ntc, "logcounts")[dx_res_sig$X, ]
sum_expression_sig <- colSums(expr_data_sig)

dx_res_fil = dx_res[which(dx_res$p_value_TRUE < 0.05),]
expr_data_fil <- assay(spe_ntc, "logcounts")[dx_res_fil$X, ]
sum_expression_fil <- colSums(expr_data_fil)

pdf(here("plots", "image_processing", "enrichment", "neuropil_Boyi.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = neuropil_pos, y = sum_expression_fil)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "Sum of Expression of -DAPI DEGs",
       x = "-DAPI",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = neuropil_pos, y = sum_expression_sig)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title = "Sum of Expression of neuropil DEGs",
       x = "-DAPI",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(logFC_TRUE ~ S_logFC, data=dx_res,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$S_padj < 0.05),],
	main = paste(n_sig, "Neuropil Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05),],
	main = paste(n_fil, "-DAPI Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),],
	main = paste(n_both, "Neuropil & -DAPI Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuropil Enrichment (log2FC)",
	ylab = "-DAPI Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_dx_res = dx_res %>% filter(p_value_TRUE < 0.05) %>% arrange(p_value_TRUE)
head(filtered_dx_res,20)

#                 X t_stat_FALSE t_stat_TRUE p_value_FALSE p_value_TRUE
#1  ENSG00000127585    -8.555610    8.555610  4.737761e-12 4.737761e-12
#2  ENSG00000078369    -8.381838    8.381838  9.416218e-12 9.416218e-12
#3  ENSG00000010404    -7.931712    7.931712  5.610576e-11 5.610576e-11
#4  ENSG00000154118    -7.427272    7.427272  4.165978e-10 4.165978e-10
#5  ENSG00000162545    -7.344170    7.344170  5.796201e-10 5.796201e-10
#6  ENSG00000167658    -7.284020    7.284020  7.360843e-10 7.360843e-10
#7  ENSG00000129682    -7.233495    7.233495  8.996608e-10 8.996608e-10
#8  ENSG00000129355    -7.224425    7.224425  9.326548e-10 9.326548e-10
#9  ENSG00000100167    -7.189033    7.189033  1.073372e-09 1.073372e-09
#10 ENSG00000272899    -7.184610    7.184610  1.092382e-09 1.092382e-09
#11 ENSG00000241360    -7.122447    7.122447  1.398060e-09 1.398060e-09
#12 ENSG00000136002    -7.121766    7.121766  1.401841e-09 1.401841e-09
#13 ENSG00000198944    -7.084892    7.084892  1.622684e-09 1.622684e-09
#14 ENSG00000078018    -7.039217    7.039217  1.944925e-09 1.944925e-09
#15 ENSG00000119487    -6.971462    6.971462  2.544096e-09 2.544096e-09
#16 ENSG00000167733    -6.959984    6.959984  2.662458e-09 2.662458e-09
#17 ENSG00000003249    -6.836876    6.836876  4.334428e-09 4.334428e-09
#18 ENSG00000185246     6.812598   -6.812598  4.771213e-09 4.771213e-09
#19 ENSG00000130956    -6.798116    6.798116  5.052364e-09 5.052364e-09
#20 ENSG00000110697    -6.792818    6.792818  5.159295e-09 5.159295e-09
#      fdr_FALSE     fdr_TRUE logFC_FALSE logFC_TRUE         ensembl      gene
#1  7.399735e-08 7.399735e-08  -0.5274917  0.5274917 ENSG00000127585    FBXL16
#2  7.399735e-08 7.399735e-08  -0.2786003  0.2786003 ENSG00000078369      GNB1
#3  2.939381e-07 2.939381e-07  -0.4126411  0.4126411 ENSG00000010404       IDS
#4  1.636917e-06 1.636917e-06  -0.4229184  0.4229184 ENSG00000154118      JPH3
#5  1.716897e-06 1.716897e-06  -0.4415287  0.4415287 ENSG00000162545   CAMK2N1
#6  1.716897e-06 1.716897e-06  -0.4004307  0.4004307 ENSG00000167658      EEF2
#7  1.716897e-06 1.716897e-06  -0.3165889  0.3165889 ENSG00000129682     FGF13
#8  1.716897e-06 1.716897e-06  -0.2590204  0.2590204 ENSG00000129355    CDKN2D
#9  1.716897e-06 1.716897e-06  -0.2939799  0.2939799 ENSG00000100167   SEPTIN3
#10 1.716897e-06 1.716897e-06  -0.5188502  0.5188502 ENSG00000272899 ATP6V1FNB
#11 1.836062e-06 1.836062e-06  -0.3657522  0.3657522 ENSG00000241360      PDXP
#12 1.836062e-06 1.836062e-06  -0.2532224  0.2532224 ENSG00000136002   ARHGEF4
#13 1.961825e-06 1.961825e-06  -0.5970234  0.5970234 ENSG00000198944    SOWAHA
#14 2.183457e-06 2.183457e-06  -0.3937867  0.3937867 ENSG00000078018      MAP2
#15 2.615366e-06 2.615366e-06  -0.2490534  0.2490534 ENSG00000119487   MAPKAP1
#16 2.615366e-06 2.615366e-06  -0.3722432  0.3722432 ENSG00000167733  HSD11B1L
#17 4.007306e-06 4.007306e-06  -0.4203311  0.4203311 ENSG00000003249    DBNDD1
#18 4.054432e-06 4.054432e-06   0.3394436 -0.3394436 ENSG00000185246    PRPF39
#19 4.054432e-06 4.054432e-06  -0.2981418  0.2981418 ENSG00000130956     HABP4
#20 4.054432e-06 4.054432e-06  -0.4113081  0.4113081 ENSG00000110697   PITPNM1



### marker genes ##
library(scater)
neuropil_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "neuropil_pseudo.rds"))
rownames(neuropil_pseudo) = rowData(neuropil_pseudo)$gene_name

p = plotExpression(neuropil_pseudo, c("CAMK2A"), x = "neuropil_pos", exprs_values = "logcounts")
ggsave(here("plots", "image_processing", "enrichment", "CAMK2A.png"), plot = p, width = 8, height = 6, dpi = 300)
p = plotExpression(neuropil_pseudo, c("DLG4"), x = "neuropil_pos", exprs_values = "logcounts")
ggsave(here("plots", "image_processing", "enrichment", "DLG4.png"), plot = p, width = 8, height = 6, dpi = 300)

### volcanoplots
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "neuropil_dx_res.csv"))
dx_res$sig <- with(dx_res, 
	                   (logFC_TRUE > log2(1.1) & fdr_TRUE < 0.02) | 
	                   (logFC_TRUE < -log2(1.1) & fdr_TRUE < 0.02)
	)
genes = read_excel(here("code/image_processing/enrichment/Maddy_SPG_volcano_highlightDEGs.xlsx"))

	dx_res <- dx_res %>%
	  mutate(label = ifelse(gene %in% genes$neuropil, gene, NA)) # Only label matching genes

    	dx_res = dx_res %>% mutate(sig = ifelse(gene %in% genes$neuropil, TRUE, FALSE))
		
    	# Create the plot
    	p <- ggplot(data = dx_res, aes(x = logFC_TRUE, y = -log10(p_value_TRUE), color = sig)) +
    	  geom_point() + 
    	  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), breaks = c("TRUE", "FALSE"), name = "FDR<0.1") + 
    	  geom_label_repel(
    	    aes(label = label), 
    	    box.padding = 0.35, 
    	    point.padding = 0.5, 
    	    segment.color = "red",  # Red lines connecting labels to points
    	    segment.size = 3,     # Thickness of the lines
    	    label.size = 0.5,         # Font size of labels
    	    color = "red",          # Label color
    	    fill = "white",         # Fill color of the bounding box
    	    fontface = "bold",      # Make the text bold
    	    max.overlaps = Inf,  # Allow all labels to be shown even if overlapping
    	    na.rm = FALSE
    	  ) +
    	  geom_hline(yintercept = -log10(0.1), color = "black", linetype = "dashed", size = 1) +
    	  theme_minimal() +
    	  labs(
    	    x = expression(log[2]~"fold change"), 
    	    y = expression(-log[10]~FDR), 
    	    title = "Differential expression in neuropil+ spots"
    	  )
			  
ggsave(here("plots", "image_processing", "enrichment", "neuropil_volcano.pdf"), plot = p, width = 8, height = 6, dpi = 300)
