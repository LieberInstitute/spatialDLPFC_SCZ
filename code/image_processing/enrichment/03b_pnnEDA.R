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
load(file = here("processed-data", "image_processing", "enrichment", "pnn_dx_res_andrew_sub.Rdata"))
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "pnn_dx_res.csv"))
pnn = read_excel(here("raw-data/images/SPG_Spot_Valid/PNN/Supple_Table_DataS4_PNN.xlsx"), sheet = "PNN Energy")
pnn = as.data.frame(pnn)

spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)

### Andrew stats 
#table(neuropil$padj < 0.05, sign(neuropil$stat))
stats$S_padj = pnn[match(stats$Symbol, toupper(pnn$gene_acronym)), 8]
stats$S_FDR = pnn[match(stats$Symbol, toupper(pnn$gene_acronym)), 7]
stats$FDR <- p.adjust(stats$P.Value, method = "BH")

ct = cor.test(stats$FDR, stats$S_FDR)
ct_sig = with(stats[which(stats$S_padj < 0.05),], cor.test(FDR, S_FDR))
ct_fil = with(stats[which(stats$P.Value < 0.05),], cor.test(FDR, S_FDR))
ct_both = with(stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),], cor.test(FDR, S_FDR))

n_all = sum(complete.cases(stats[,c("FDR", "S_FDR")]))
n_sig = sum(complete.cases(stats[which(stats$S_padj < 0.05),c("FDR", "S_FDR")]))
n_fil = sum(complete.cases(stats[which(stats$P.Value < 0.05),c("FDR", "S_FDR")]))
n_both = sum(complete.cases(stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),c("FDR", "S_FDR")]))

stats_sig = stats[which(stats$P.Value < 0.05),]
expr_data <- assay(spe_ntc, "logcounts")[rownames(stats_sig), ]
sum_expression <- colSums(expr_data)

pdf(here("plots", "image_processing", "enrichment", "pnn_andrewSub.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = pnn_pos, y = sum_expression)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA)
  labs(title = "Sum of Expression of wfa DEGs",
       x = "wfa",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(FDR ~ S_FDR, data=stats,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(FDR ~ S_FDR, 
	data=stats[which(stats$S_padj < 0.05),],
	main = paste(n_sig, "pnn Significant genes"),
	pch = 21, bg = "grey",
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(FDR ~ S_FDR, 
	data=stats[which(stats$P.Value < 0.05),],
	main = paste(n_fil, "wfa Significant genes"),
	pch = 21, bg = "grey",
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(FDR ~ S_FDR, 
	data=stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),],
	main = paste(n_both, "pnn & wfa Significant genes"),
	pch = 21, bg = "grey",
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_stats = stats %>% filter(P.Value < 0.05) %>% arrange(P.Value)
head(filtered_stats,20)

#                      logFC         t       P.Value   Symbol       S_padj
#ENSG00000100362  0.27634102  23.19024 5.033138e-117    PVALB 5.357403e-11
#ENSG00000100285  0.37354366  22.15484 3.940991e-107     NEFH 4.986080e+00
#ENSG00000101439 -0.41608316 -16.82209  5.910553e-63     CST3 1.591816e+02
#ENSG00000110436 -0.51673170 -16.51260  9.575404e-61   SLC1A2 3.261878e+02
#ENSG00000006128  0.18513537  16.41081  5.003423e-60     TAC1           NA
#ENSG00000104722  0.38396332  14.77404  4.567038e-49     NEFM 5.795977e-10
#ENSG00000139190  0.16799317  14.66376  2.281793e-48    VAMP1 8.633578e-16
#ENSG00000143858  0.07242854  14.59231  6.430547e-48     SYT2 8.274130e-13
#ENSG00000087250 -0.24986330 -13.68307  2.231659e-42      MT3 1.704377e+02
#ENSG00000135821 -0.36525405 -13.19353  1.544934e-39     GLUL 3.416513e-12
#ENSG00000143153  0.29620691  13.15114  2.692372e-39   ATP1B1 2.165318e+03
#ENSG00000131711  0.28226613  12.98386  2.370345e-38    MAP1B 1.612398e+03
#ENSG00000277586  0.27978304  12.05082  2.685195e-33     NEFL 1.852697e-14
#ENSG00000130203 -0.30651063 -11.94841  9.151244e-33     APOE 3.697256e-02
#ENSG00000173267  0.22321201  11.75836  8.670598e-32     SNCG 6.270298e-02
#ENSG00000079215 -0.27883095 -11.07641  2.075889e-28   SLC1A3 7.635966e-09
#ENSG00000014641  0.25663921  11.02463  3.679705e-28     MDH1 8.543665e+02
#ENSG00000080824  0.22737847  10.95683  7.755658e-28 HSP90AA1 8.381050e+02
#ENSG00000152595  0.05250463  10.61316  3.170336e-26     MEPE 5.890478e+03
#ENSG00000018625 -0.26586627 -10.54212  6.729858e-26   ATP1A2 1.174156e+02

### Boyi stats 
dx_res$S_FDR = pnn[match(dx_res$gene, toupper(pnn$gene_acronym)), 7]
dx_res$S_padj = pnn[match(dx_res$gene, toupper(pnn$gene_acronym)), 8]

ct = cor.test(dx_res$fdr_TRUE, dx_res$S_FDR)
ct_sig = with(dx_res[which(dx_res$S_padj < 0.05),], cor.test(fdr_TRUE, S_FDR))
ct_fil = with(dx_res[which(dx_res$p_value_TRUE < 0.05),], cor.test(fdr_TRUE, S_FDR))
ct_both = with(dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),], cor.test(fdr_TRUE, S_FDR))

n_all = sum(complete.cases(dx_res[,c("fdr_TRUE", "S_FDR")]))
n_sig = sum(complete.cases(dx_res[which(dx_res$S_padj < 0.05),c("fdr_TRUE", "S_FDR")]))
n_fil = sum(complete.cases(dx_res[which(dx_res$p_value_TRUE < 0.05),c("fdr_TRUE", "S_FDR")]))
n_both = sum(complete.cases(dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),c("fdr_TRUE", "S_FDR")]))

dx_res_sig = dx_res[which(dx_res$S_padj < 0.05),]
expr_data_sig <- assay(spe_ntc, "logcounts")[dx_res_sig$X, ]
sum_expression_sig <- colSums(expr_data)

dx_res_fil = dx_res[which(dx_res$p_value_TRUE < 0.05),]
expr_data_fil <- assay(spe_ntc, "logcounts")[dx_res_sig$X, ]
sum_expression_fil <- colSums(expr_data)

pdf(here("plots", "image_processing", "enrichment", "pnn_Boyi.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = pnn_pos, y = sum_expression_fil)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "Sum of Expression of wfa DEGs",
       x = "wfa",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = pnn_pos, y = sum_expression_sig)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title = "Sum of Expression of pnn DEGs",
       x = "wfa",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(fdr_TRUE~ S_FDR, data=dx_res,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(fdr_TRUE~ S_FDR, 
	data=dx_res[which(dx_res$S_padj < 0.05),],
	main = paste(n_sig, "pnn Significant genes"),
	pch = 21, bg = "grey",
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(fdr_TRUE~ S_FDR, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05),],
	main = paste(n_fil, "wfa Significant genes"),
	pch = 21, bg = "grey",
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(fdr_TRUE~ S_FDR, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),],
	main = paste(n_both, "wfa & pnn Significant genes"),
	pch = 21, bg = "grey",
	xlab = "pnn Enrichment (FDR)",
	ylab = "wfa Enrichment (FDR)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_dx_res = dx_res %>% filter(p_value_TRUE < 0.05) %>% arrange(p_value_TRUE)
head(filtered_dx_res,20)
#                 X t_stat_FALSE t_stat_TRUE p_value_FALSE p_value_TRUE
#1  ENSG00000100362    -17.40599    17.40599  4.224910e-25 4.224910e-25
#2  ENSG00000153012    -15.08734    15.08734  4.597851e-22 4.597851e-22
#3  ENSG00000196482    -13.79174    13.79174  3.044077e-20 3.044077e-20
#4  ENSG00000113140     12.85849   -12.85849  7.121865e-19 7.121865e-19
#5  ENSG00000170745    -12.79631    12.79631  8.821446e-19 8.821446e-19
#6  ENSG00000029534    -12.43856    12.43856  3.051324e-18 3.051324e-18
#7  ENSG00000041982     12.39901   -12.39901  3.503503e-18 3.503503e-18
#8  ENSG00000006128    -12.27643    12.27643  5.383564e-18 5.383564e-18
#9  ENSG00000100285    -12.14415    12.14415  8.576852e-18 8.576852e-18
#10 ENSG00000145423     12.13458   -12.13458  8.871598e-18 8.871598e-18
#11 ENSG00000006611     12.01198   -12.01198  1.368921e-17 1.368921e-17
#12 ENSG00000136750    -11.96937    11.96937  1.592355e-17 1.592355e-17
#13 ENSG00000147145     11.66454   -11.66454  4.728086e-17 4.728086e-17
#14 ENSG00000182836    -11.49932    11.49932  8.569483e-17 8.569483e-17
#15 ENSG00000151136    -11.42094    11.42094  1.137565e-16 1.137565e-16
#16 ENSG00000081138    -11.41095    11.41095  1.179474e-16 1.179474e-16
#17 ENSG00000115756     11.41069   -11.41069  1.180565e-16 1.180565e-16
#18 ENSG00000166006    -11.18499    11.18499  2.681102e-16 2.681102e-16
#19 ENSG00000145147    -10.77652    10.77652  1.201361e-15 1.201361e-15
#20 ENSG00000198944     10.65849   -10.65849  1.859827e-15 1.859827e-15
#      fdr_FALSE     fdr_TRUE logFC_FALSE logFC_TRUE         ensembl   gene
#1  6.294694e-21 6.294694e-21  -2.4162749  2.4162749 ENSG00000100362  PVALB
#2  3.425169e-18 3.425169e-18  -1.7261176  1.7261176 ENSG00000153012   LGI2
#3  1.511790e-16 1.511790e-16  -1.1909392  1.1909392 ENSG00000196482  ESRRG
#4  2.628614e-15 2.628614e-15   1.5048665 -1.5048665 ENSG00000113140  SPARC
#5  2.628614e-15 2.628614e-15  -1.7788742  1.7788742 ENSG00000170745  KCNS3
#6  7.456956e-15 7.456956e-15  -1.3116708  1.3116708 ENSG00000029534   ANK1
#7  7.456956e-15 7.456956e-15   3.9437456 -3.9437456 ENSG00000041982    TNC
#8  1.002622e-14 1.002622e-14  -1.7830258  1.7830258 ENSG00000006128   TAC1
#9  1.321779e-14 1.321779e-14  -1.7699511  1.7699511 ENSG00000100285   NEFH
#10 1.321779e-14 1.321779e-14   4.1214362 -4.1214362 ENSG00000145423  SFRP2
#11 1.854141e-14 1.854141e-14   4.2973722 -4.2973722 ENSG00000006611  USH1C
#12 1.977041e-14 1.977041e-14  -1.2353117  1.2353117 ENSG00000136750   GAD2
#13 5.418751e-14 5.418751e-14   3.7151227 -3.7151227 ENSG00000147145  LPAR4
#14 9.119766e-14 9.119766e-14  -0.9995345  0.9995345 ENSG00000182836 PLCXD3
#15 1.034661e-13 1.034661e-13  -1.3612362  1.3612362 ENSG00000151136 BTBD11
#16 1.034661e-13 1.034661e-13  -1.2895433  1.2895433 ENSG00000081138   CDH7
#17 1.034661e-13 1.034661e-13   1.3489496 -1.3489496 ENSG00000115756 HPCAL1
#18 2.219207e-13 2.219207e-13  -1.0025162  1.0025162 ENSG00000166006  KCNC2
#19 9.420565e-13 9.420565e-13  -0.8714369  0.8714369 ENSG00000145147  SLIT2
#20 1.385478e-12 1.385478e-12   0.8902372 -0.8902372 ENSG00000198944 SOWAHA


### marker genes

library(scater)
pnn_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "pnn_pseudo.rds"))
rownames(pnn_pseudo) = rowData(pnn_pseudo)$gene_name

p = plotExpression(pnn_pseudo, c("PVALB"), x = "pnn_pos", exprs_values = "logcounts")

ggsave(here("plots", "image_processing", "enrichment", "pnn_genes.png"), plot = p, width = 8, height = 6, dpi = 300)

### volcanoplots
library(ggplot2)
library(grid) 
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "pnn_dx_res.csv"))
#dx_res$sig <- with(dx_res, 
#	                   (logFC_TRUE > log2(1.2) & fdr_TRUE < 0.1) | 
#	                   (logFC_TRUE < -log2(1.2) & fdr_TRUE < 0.1)
#	)
genes = read_excel(here("code/image_processing/enrichment/Maddy_SPG_volcano_highlightDEGs.xlsx"))
genes$wfa[genes$wfa == "HPLN4"] = "HAPLN4"
genes$wfa[genes$wfa == "HPLN1"] = "HAPLN1"
	dx_res <- dx_res %>%
	  mutate(label = ifelse(gene %in% genes$wfa, gene, NA)) # Only label matching genes
	
	dx_res = dx_res %>% mutate(sig = ifelse(gene %in% genes$wfa, TRUE, FALSE))
		
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
	    title = "Differential expression in wfa+ spots"
	  )
			  
ggsave(here("plots", "image_processing", "enrichment", "pnn_volcano.pdf"), plot = p, width = 8, height = 6, dpi = 300)

