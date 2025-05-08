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
load(file = here("processed-data", "image_processing", "enrichment", "neun_dx_res_andrew_sub.Rdata"))
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "neun_dx_res.csv"))
neun1 = read_excel(here("raw-data/images/SPG_Spot_Valid/NeuN/TableS13_marker_stats_supp_table.xlsx"), sheet="marker_stats_supp_table")
neun1 = as.data.frame(neun1)
neun1$padj = as.numeric(neun1$log.p.value)
celltypes = c("Inhib", "Excit", "Excit_L3/4/5", "Excit_L3", "Excit_L4", "Excit_L6", "Excit_L5/6", "Excit_L5", "Excit_L2/3")
neun = neun1[neun1$cellType.target %in% celltypes,]

spe_ntc = readRDS(here("processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
colData(spe_ntc)$slide_id <- sapply(strsplit(colData(spe_ntc)$sample_id, "_"), `[`, 1)

### Andrew stats 
#table(neuropil$padj < 0.05, sign(neuropil$stat))
stats$S_logFC = neun$std.logFC[match(rownames(stats), neun$gene)]
stats$S_padj = neun$log.p.value[match(rownames(stats), neun$gene)]

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

pdf(here("plots", "image_processing", "enrichment", "neun_andrewSub.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = neun_pos, y = sum_expression)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA)
  labs(title = "Sum of Expression of neun DEGs",
       x = "NeuN",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(logFC ~ S_logFC, data=stats,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$S_padj < 0.05),],
	main = paste(n_sig, "Neuron Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05),],
	main = paste(n_fil, "NeuN Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC ~ S_logFC, 
	data=stats[which(stats$P.Value < 0.05 & stats$S_padj < 0.05),],
	main = paste(n_both, "NeuN & Vasc Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_stats = stats %>% filter(P.Value < 0.05) %>% arrange(P.Value)
head(filtered_stats,20)

#                     logFC         t       P.Value    Symbol  S_logFC    S_padj
#ENSG00000198888 -0.4514841 -39.91653  0.000000e+00    MT-ND1       NA        NA
#ENSG00000198763 -0.4663654 -42.62854  0.000000e+00    MT-ND2       NA        NA
#ENSG00000198804 -0.3426674 -40.58548  0.000000e+00    MT-CO1       NA        NA
#ENSG00000198712 -0.4307839 -40.49725  0.000000e+00    MT-CO2       NA        NA
#ENSG00000198899 -0.4101196 -40.03602  0.000000e+00   MT-ATP6       NA        NA
#ENSG00000198938 -0.3829693 -39.67242  0.000000e+00    MT-CO3       NA        NA
#ENSG00000198840 -0.3979028 -41.55608  0.000000e+00    MT-ND3       NA        NA
#ENSG00000198886 -0.3731548 -38.85530 1.539376e-315    MT-ND4       NA        NA
#ENSG00000269028 -0.4985691 -37.53056 2.098538e-295 MTRNR2L12 1.011274 -1385.197
#ENSG00000198727 -0.4039225 -37.46769 1.832509e-294    MT-CYB       NA        NA
#ENSG00000080824  0.4457340  32.99542 5.708059e-231  HSP90AA1 0.863950 -1019.018
#ENSG00000198786 -0.4315559 -32.91642 6.542207e-230    MT-ND5       NA        NA
#ENSG00000067715  0.4557444  29.06237 5.789788e-181      SYT1       NA        NA
#ENSG00000014641  0.3801338  24.72824 1.704108e-132      MDH1       NA        NA
#ENSG00000131711  0.3507539  24.38379 5.971022e-129     MAP1B 1.697614 -3728.574
#ENSG00000104722  0.4180418  24.26812 9.038535e-128      NEFM       NA        NA
#ENSG00000228253 -0.3606274 -23.21901 2.635868e-117   MT-ATP8       NA        NA
#ENSG00000212907 -0.3360731 -22.95113 1.057051e-114   MT-ND4L       NA        NA
#ENSG00000164830  0.3351364  22.87101 6.268111e-114      OXR1       NA        NA
#ENSG00000139970  0.3386155  22.61573 1.752458e-111      RTN1       NA        NA

### Boyi stats 
dx_res$S_logFC = neun$std.logFC[match(dx_res$X, neun$gene)]
dx_res$S_padj = neun$log.p.value[match(dx_res$X, neun$gene)]

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
sum_expression_sig <- colSums(expr_data)

dx_res_fil = dx_res[which(dx_res$p_value_TRUE < 0.05),]
expr_data_fil <- assay(spe_ntc, "logcounts")[dx_res_sig$X, ]
sum_expression_fil <- colSums(expr_data)

pdf(here("plots", "image_processing", "enrichment", "neun_Boyi.pdf"))
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = neun_pos, y = sum_expression_fil)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "Sum of Expression of neun DEGs",
       x = "NeuN",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
p = ggplot(as.data.frame(colData(spe_ntc)), aes(x = neun_pos, y = sum_expression_sig)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title = "Sum of Expression of neuron DEGs",
       x = "NeuN",
       y = "Sum of Expression") +
  theme_minimal()
print(p)
plot(logFC_TRUE ~ S_logFC, data=dx_res,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed genes"),
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$S_padj < 0.05),],
	main = paste(n_sig, "Neuron Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05),],
	main = paste(n_fil, "NeuN Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_fil$estimate, 3)),cex=1.5)
plot(logFC_TRUE ~ S_logFC, 
	data=dx_res[which(dx_res$p_value_TRUE < 0.05 & dx_res$S_padj < 0.05),],
	main = paste(n_both, "NeuN & Vasc Significant genes"),
	pch = 21, bg = "grey",
	xlab = "Neuron Enrichment (log2FC)",
	ylab = "NeuN Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct_both$estimate, 3)),cex=1.5)
dev.off()

filtered_dx_res = dx_res %>% filter(p_value_TRUE < 0.05) %>% arrange(p_value_TRUE)
head(filtered_dx_res,20)
#                 X t_stat_FALSE t_stat_TRUE p_value_FALSE p_value_TRUE
#1  ENSG00000113140     17.91761   -17.91761  5.207234e-26 5.207234e-26
#2  ENSG00000249306     14.74467   -14.74467  8.514256e-22 8.514256e-22
#3  ENSG00000120693     14.11863   -14.11863  6.675217e-21 6.675217e-21
#4  ENSG00000251442     13.81273   -13.81273  1.859139e-20 1.859139e-20
#5  ENSG00000119139     13.75690   -13.75690  2.244161e-20 2.244161e-20
#6  ENSG00000124942     12.91323   -12.91323  4.050873e-19 4.050873e-19
#7  ENSG00000152558     12.91122   -12.91122  4.079294e-19 4.079294e-19
#8  ENSG00000134817     12.89842   -12.89842  4.265356e-19 4.265356e-19
#9  ENSG00000144579     12.77089   -12.77089  6.659457e-19 6.659457e-19
#10 ENSG00000130303     12.77045   -12.77045  6.669755e-19 6.669755e-19
#11 ENSG00000204592     12.65953   -12.65953  9.843006e-19 9.843006e-19
#12 ENSG00000173193     12.58255   -12.58255  1.290735e-18 1.290735e-18
#13 ENSG00000150093     12.56785   -12.56785  1.359426e-18 1.359426e-18
#14 ENSG00000213398     12.54180   -12.54180  1.490305e-18 1.490305e-18
#15 ENSG00000173905     12.48341   -12.48341  1.831934e-18 1.831934e-18
#16 ENSG00000117115     12.41710   -12.41710  2.317006e-18 2.317006e-18
#17 ENSG00000137959     12.13485   -12.13485  6.336755e-18 6.336755e-18
#18 ENSG00000204580     12.07670   -12.07670  7.806013e-18 7.806013e-18
#19 ENSG00000177469     11.97336   -11.97336  1.131935e-17 1.131935e-17
#20 ENSG00000147145     11.87244   -11.87244  1.629315e-17 1.629315e-17
#      fdr_FALSE     fdr_TRUE logFC_FALSE logFC_TRUE         ensembl      gene
#1  8.142553e-22 8.142553e-22   1.5825707 -1.5825707 ENSG00000113140     SPARC
#2  6.656871e-18 6.656871e-18   2.1269116 -2.1269116 ENSG00000249306 LINC01411
#3  3.479345e-17 3.479345e-17   1.6827160 -1.6827160 ENSG00000120693     SMAD9
#4  7.018389e-17 7.018389e-17   1.6097174 -1.6097174 ENSG00000251442 LINC01094
#5  7.018389e-17 7.018389e-17   1.8343157 -1.8343157 ENSG00000119139      TJP2
#6  8.337172e-16 8.337172e-16   1.5153495 -1.5153495 ENSG00000124942     AHNAK
#7  8.337172e-16 8.337172e-16   1.2183608 -1.2183608 ENSG00000152558   TMEM123
#8  8.337172e-16 8.337172e-16   3.5599062 -3.5599062 ENSG00000134817     APLNR
#9  1.042950e-15 1.042950e-15   1.0631375 -1.0631375 ENSG00000144579    CTDSP1
#10 1.042950e-15 1.042950e-15   1.3521801 -1.3521801 ENSG00000130303      BST2
#11 1.399228e-15 1.399228e-15   1.2142711 -1.2142711 ENSG00000204592     HLA-E
#12 1.635180e-15 1.635180e-15   1.1101392 -1.1101392 ENSG00000173193    PARP14
#13 1.635180e-15 1.635180e-15   0.9097106 -0.9097106 ENSG00000150093     ITGB1
#14 1.664565e-15 1.664565e-15   1.7373679 -1.7373679 ENSG00000213398      LCAT
#15 1.909730e-15 1.909730e-15   1.0105647 -1.0105647 ENSG00000173905    GOLIM4
#16 2.264439e-15 2.264439e-15   1.3255514 -1.3255514 ENSG00000117115     PADI2
#17 5.828696e-15 5.828696e-15   1.5336329 -1.5336329 ENSG00000137959    IFI44L
#18 6.781257e-15 6.781257e-15   1.2580577 -1.2580577 ENSG00000204580      DDR1
#19 9.315822e-15 9.315822e-15   1.6339418 -1.6339418 ENSG00000177469    CAVIN1
#20 1.273880e-14 1.273880e-14   2.4898095 -2.4898095 ENSG00000147145     LPAR4

### marker genes

library(scater)
neun_pseudo = readRDS(here("processed-data", "image_processing", "enrichment", "neun_pseudo.rds"))
rownames(neun_pseudo) = rowData(neun_pseudo)$gene_name

p = plotExpression(neun_pseudo, c("SNAP25"), x = "neun_pos", exprs_values = "logcounts")
ggsave(here("plots", "image_processing", "enrichment", "SNAP25.png"), plot = p, width = 8, height = 6, dpi = 300)
p = plotExpression(neun_pseudo, c("RBFOX3"), x = "neun_pos", exprs_values = "logcounts")
ggsave(here("plots", "image_processing", "enrichment", "RBFOX3.png"), plot = p, width = 8, height = 6, dpi = 300)


### volcanoplots
dx_res = read.csv(file = here("processed-data", "image_processing", "enrichment", "neun_dx_res.csv"))
dx_res$sig <- with(dx_res, 
	                   (logFC_TRUE > log2(1.2) & fdr_TRUE < 0.0001) | 
	                   (logFC_TRUE < -log2(1.2) & fdr_TRUE < 0.0001)
	)
genes = read_excel(here("code/image_processing/enrichment/Maddy_SPG_volcano_highlightDEGs.xlsx"))

	dx_res <- dx_res %>%
	  mutate(label = ifelse(gene %in% genes$neun, gene, NA)) # Only label matching genes

  	dx_res = dx_res %>% mutate(sig = ifelse(gene %in% genes$neun, TRUE, FALSE))
		
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
  	    title = "Differential expression in neun+ spots"
  	  )
			  
ggsave(here("plots", "image_processing", "enrichment", "neun_volcano.png"), plot = p, width = 8, height = 6, dpi = 300)

