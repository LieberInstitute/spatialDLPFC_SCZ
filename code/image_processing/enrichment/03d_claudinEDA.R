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
  library(ggrepel)
})

# Load Data ----
load(file = here("processed-data", "image_processing", "enrichment", "vasc_dx_res_andrew_sub.Rdata"))
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "vasc_dx_res.csv"))
vasc = read_excel(here("raw-data/images/SPG_Spot_Valid/Vasculature/Supple_Table_2_Vas.xlsx"), sheet="Post Mortem Vascular Subcluster")
vasc = as.data.frame(vasc)
vasc$padj = as.numeric(vasc$p_val_adj)
spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)

### Andrew stats 
#table(neuropil$padj < 0.05, sign(neuropil$stat))
stats$S_logFC = vasc$avg_log2FC[match(stats$Symbol, vasc$Column1)]
stats$S_padj = vasc$p_val_adj[match(stats$Symbol, vasc$Column1)]

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

pdf(here("plots", "image_processing", "enrichment", "vasc_andrewSub.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = vasc_pos, y = sum_expression)) +
  geom_violin() + geom_boxplot(width = 0.1) +
  labs(title = "Sum of Expression of claudin5 DEGs",
       x = "Claudin5",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(logFC ~ S_logFC, data=stats,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$S_padj < 0.05),],
	main = paste(n_sig, "Vasculature Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05),],
	main = paste(n_fil, "Claudin5 Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),],
	main = paste(n_both, "Claudin5 & Vasc Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_stats = stats %>% filter(P.Value < 0.05) %>% arrange(P.Value)
head(filtered_stats,20)
#                      logFC         t      P.Value  Symbol       S_logFC S_padj
#ENSG00000101335  0.25595323  18.96472 2.520141e-79    MYL9  0.000000e+00  0.184
#ENSG00000175899  0.41922461  18.92471 5.298743e-79     A2M  2.153690e-16  0.428
#ENSG00000137831  0.27168234  18.53561 6.752975e-76    UACA 2.868688e-207  0.663
#ENSG00000122786  0.36692032  17.65971 3.934254e-69   CALD1  1.043255e-82  0.679
#ENSG00000054654  0.21529625  16.71786 3.312364e-62   SYNE2  3.259449e-10  0.702
#ENSG00000149591  0.23633226  16.51968 8.532437e-61   TAGLN  0.000000e+00  0.135
#ENSG00000107796  0.18054435  16.35939 1.149334e-59   ACTA2  0.000000e+00  0.162
#ENSG00000198899 -0.33797991 -15.86115 3.186566e-56 MT-ATP6  1.298004e-43  0.833
#ENSG00000184113  0.24791886  15.75917 1.568322e-55   CLDN5  2.698333e-40  0.355
#ENSG00000198840 -0.31180873 -15.59548 1.983105e-54  MT-ND3  4.874818e-48  0.696
#ENSG00000177469  0.17130651  15.31950 1.349193e-52  CAVIN1  2.623272e-09  0.197
#ENSG00000198763 -0.34630451 -15.11377 2.991011e-51  MT-ND2  7.214917e-53  0.578
#ENSG00000133392  0.16921581  14.90048 7.119504e-50   MYH11  0.000000e+00  0.030
#ENSG00000111341  0.13243166  14.58072 7.603073e-48     MGP 9.358376e-122  0.123
#ENSG00000198888 -0.34169884 -14.50846 2.155417e-47  MT-ND1  2.920430e-57  0.672
#ENSG00000198467  0.15578070  14.46121 4.248918e-47    TPM2 1.090784e-279  0.109
#ENSG00000169908  0.11810911  14.40200 9.914974e-47  TM4SF1  1.604119e-45  0.212
#ENSG00000168497  0.16486067  14.35138 2.040870e-46  CAVIN2  0.000000e+00  0.105
#ENSG00000196154  0.09490481  14.25152 8.416020e-46  S100A4  0.000000e+00  0.015
#ENSG00000198938 -0.28366428 -14.11519 5.733471e-45  MT-CO3  1.549489e-68  0.855

### Boyi stats 
dx_res$S_logFC = vasc$avg_log2FC[match(dx_res$gene, vasc$Column1)]
dx_res$S_padj = vasc$p_val_adj[match(dx_res$gene, vasc$Column1)]

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
expr_data_fil <- assay(spe_ntc, "logcounts")[dx_res_sig$X, ]
sum_expression_fil <- colSums(expr_data_sig)


pdf(here("plots", "image_processing", "enrichment", "vasc_Boyi.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = vasc_pos, y = sum_expression_fil)) +
  geom_violin() + geom_boxplot(width = 0.1) +
  labs(title = "Sum of Expression of claudin5 DEGs",
       x = "Claudin5",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = vasc_pos, y = sum_expression_sig)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "Sum of Expression of vasc DEGs",
       x = "Claudin5",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(logFC_TRUE ~ S_logFC, data=dx_res,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$S_padj < 0.05),],
	main = paste(n_sig, "Vasculature Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05),],
	main = paste(n_fil, "Claudin5 Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),],
	main = paste(n_both, "Claudin5 & Vasc Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Vasculature Enrichment (log2FC)",
	ylab = "Claudin5 Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_dx_res = dx_res %>% filter(p_value_TRUE < 0.05) %>% arrange(p_value_TRUE)
head(filtered_dx_res,20)
#                 X t_stat_FALSE t_stat_TRUE p_value_FALSE p_value_TRUE
#1  ENSG00000123243    -17.03342    17.03342  1.146660e-24 1.146660e-24
#2  ENSG00000163072    -15.17277    15.17277  3.277544e-22 3.277544e-22
#3  ENSG00000138795    -14.55596    14.55596  2.344677e-21 2.344677e-21
#4  ENSG00000054654    -14.48838    14.48838  2.917185e-21 2.917185e-21
#5  ENSG00000069122    -14.34785    14.34785  4.603321e-21 4.603321e-21
#6  ENSG00000118777    -13.94841    13.94841  1.706447e-20 1.706447e-20
#7  ENSG00000162618    -13.73236    13.73236  3.495782e-20 3.495782e-20
#8  ENSG00000175899    -13.62321    13.62321  5.033601e-20 5.033601e-20
#9  ENSG00000143248    -13.53724    13.53724  6.715324e-20 6.715324e-20
#10 ENSG00000127329    -13.23197    13.23197  1.883110e-19 1.883110e-19
#11 ENSG00000213949    -13.12804    13.12804  2.682363e-19 2.682363e-19
#12 ENSG00000102755    -13.00499    13.00499  4.085069e-19 4.085069e-19
#13 ENSG00000085563    -12.91389    12.91389  5.584819e-19 5.584819e-19
#14 ENSG00000078596    -12.90409    12.90409  5.776130e-19 5.776130e-19
#15 ENSG00000168497    -12.70142    12.70142  1.162850e-18 1.162850e-18
#16 ENSG00000082438    -12.53367    12.53367  2.083325e-18 2.083325e-18
#17 ENSG00000164035    -12.40772    12.40772  3.235369e-18 3.235369e-18
#18 ENSG00000163513    -12.26795    12.26795  5.285687e-18 5.285687e-18
#19 ENSG00000267107    -12.21904    12.21904  6.279798e-18 6.279798e-18
#20 ENSG00000054598    -12.17933    12.17933  7.224598e-18 7.224598e-18
#      fdr_FALSE     fdr_TRUE logFC_FALSE logFC_TRUE         ensembl    gene
#1  1.691897e-20 1.691897e-20   -1.856286   1.856286 ENSG00000123243   ITIH5
#2  2.418008e-18 2.418008e-18   -1.853839   1.853839 ENSG00000163072 NOSTRIN
#3  1.076077e-17 1.076077e-17   -1.875807   1.875807 ENSG00000138795    LEF1
#4  1.076077e-17 1.076077e-17   -1.716562   1.716562 ENSG00000054654   SYNE2
#5  1.358440e-17 1.358440e-17   -1.518678   1.518678 ENSG00000069122  ADGRF5
#6  4.196438e-17 4.196438e-17   -1.717208   1.717208 ENSG00000118777   ABCG2
#7  7.368610e-17 7.368610e-17   -1.760579   1.760579 ENSG00000162618  ADGRL4
#8  9.283847e-17 9.283847e-17   -1.629899   1.629899 ENSG00000175899     A2M
#9  1.100940e-16 1.100940e-16   -1.244594   1.244594 ENSG00000143248    RGS5
#10 2.778528e-16 2.778528e-16   -1.369440   1.369440 ENSG00000127329   PTPRB
#11 3.598025e-16 3.598025e-16   -1.742280   1.742280 ENSG00000213949   ITGA1
#12 5.022933e-16 5.022933e-16   -1.810324   1.810324 ENSG00000102755    FLT1
#13 6.087628e-16 6.087628e-16   -1.795887   1.795887 ENSG00000085563   ABCB1
#14 6.087628e-16 6.087628e-16   -1.266474   1.266474 ENSG00000078596   ITM2A
#15 1.143857e-15 1.143857e-15   -1.875659   1.875659 ENSG00000168497  CAVIN2
#16 1.921216e-15 1.921216e-15   -1.443708   1.443708 ENSG00000082438  COBLL1
#17 2.808110e-15 2.808110e-15   -1.714466   1.714466 ENSG00000164035    EMCN
#18 4.332795e-15 4.332795e-15   -1.618001   1.618001 ENSG00000163513  TGFBR2
#19 4.876759e-15 4.876759e-15   -2.049856   2.049856 ENSG00000267107  PCAT19
#20 5.329947e-15 5.329947e-15   -1.971997   1.971997 ENSG00000054598   FOXC1


### marker genes

library(scater)
vasc_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "vasc_pseudo.rds"))
rownames(vasc_pseudo) = rowData(vasc_pseudo)$gene_name

p = plotExpression(vasc_pseudo, c("CLDN5"), x = "vasc_pos", exprs_values = "logcounts")
ggsave(here("plots", "image_processing", "enrichment", "CLDN5.png"), plot = p, width = 8, height = 6, dpi = 300)


## valcano plots
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "vasc_dx_res.csv"))
dx_res$sig <- with(dx_res, 
	                   (logFC_TRUE > log2(2) & fdr_TRUE < 0.05) | 
	                   (logFC_TRUE < -log2(2) & fdr_TRUE < 0.05)
	)
genes = read_excel(here("code/image_processing/enrichment/Maddy_SPG_volcano_highlightDEGs.xlsx"))

	dx_res <- dx_res %>%
	  mutate(label = ifelse(gene %in% genes$claudin5, gene, NA)) # Only label matching genes

  	dx_res = dx_res %>% mutate(sig = ifelse(gene %in% genes$claudin5, TRUE, FALSE))
	
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
  	    na.rm = FALSE, size = 5
  	  ) +
  	  geom_hline(yintercept = -log10(0.1), color = "black", linetype = "dashed", size = 1) +
  	  theme_minimal(base_size = 18) +
  	  labs(
  	    x = expression(log[2]~"fold change"), 
  	    y = expression(-log[10]~FDR), 
  	    title = "Differential expression in Claudin5+ spots"
  	  )
			  
ggsave(here("plots", "image_processing", "enrichment", "vasc_volcano.pdf"), plot = p, width = 8, height = 6, dpi = 300)
				
