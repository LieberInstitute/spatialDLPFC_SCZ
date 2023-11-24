# library("here")
# library("ggpubr")
# library("sessioninfo")
# 
# ## Rscript -e 'sgejobs::job_single("02_explore_metrics_spaceranger", create_shell = TRUE, queue = "bluejay")'
# 
# ## Output dirs
# dir_rdata <-
#     here("processed-data", "spaceranger", "03_merged_metrics")
# dir_plots <- here("plots", "spaceranger", "03_merged_metrics")
# dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
# dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
# 
# ## Read in merged metrics
# merged_metrics <- read.csv(file.path(dir_rdata, "merged_metrics.csv"),
#     check.names = FALSE
# )
# 
# ## Compare numerical metrics one at a time across studies
# vars_class <- sapply(merged_metrics, class)
# vars_number <-
#     names(vars_class)[vars_class %in% c("integer", "numeric")]
# ## Check for change in total number of number variables
# # stopifnot(length(vars_number) == 32)
# 
# ## Define colors
# sample_status_uniq <- unique(merged_metrics$sample_status)
# # paletteer::paletteer_d("ggprism::inferno", n = 6, type = "discrete")
# # <colors>
# # 000004FF #420A68FF #932667FF #DD513AFF #FCA50AFF #FCFFA4FF
# stopifnot(length(sample_status_uniq) == 4)
# sample_status_colors <- c("#420A68FF", "#932667FF", "#DD513AFF", "#FCA50AFF")
# names(sample_status_colors) <- sample_status_uniq
# 
# ## Shorten study_name_sheet
# merged_metrics$study_name_sheet_short <-
#     gsub(
#         "spatial_|spatial|Visium_|Visium",
#         "",
#         merged_metrics$study_name_sheet
#     )
# 
# ## Use the "Sample.ID" for the "HumanPilot" samples
# i_human_pilot <- which(merged_metrics$study_name == "HumanPilot")
# merged_metrics$slide_serial_capture_area[i_human_pilot] <- merged_metrics$Sample.ID[i_human_pilot]
# 
# ## Subset to samples with metrics observed in both pre-sequencing and sequencing
# merged_metrics_both <- subset(merged_metrics, !is.na(summary_file) & !is.na(sheet_file))
# 
# ## Adapted from https://github.com/LieberInstitute/Visium_IF_AD/blob/036224121d8d56683d771cbc332ca3a40e59f004/code/03_analysis/02_sample_metrics.R#L196-L208
# ## Used https://rpkgs.datanovia.com/ggpubr/reference/ggboxplot.html
# pdf(
#     file.path(dir_plots, "spaceranger_cross_study_metrics_boxplots.pdf"),
#     useDingbats = FALSE,
#     width = length(unique(merged_metrics$study_name_sheet)) * 2.5,
#     height = 10
# )
# for (i in vars_number) {
#     # for (i in vars_number[c(1, 10, 32)]) {
#     set.seed(20220913)
#     p <-
#         ggboxplot(
#             merged_metrics[!is.na(merged_metrics[[i]]), ],
#             x = "study_name_sheet_short",
#             y = i,
#             color = "sample_status",
#             palette = sample_status_colors,
#             add = "jitter",
#             label = "slide_serial_capture_area",
#             repel = TRUE,
#             font.label = list(size = 10),
#             add.params = list(size = 3),
#             legend = "top",
#             ggtheme = theme_pubr(base_size = 20),
#             xlab = ""
#         )
#     print(p)
# }
# dev.off()
# 
# 
# ## Identify all pre-sequencing and post-sequencing numerical variables
# split_point <- which(vars_number == "Number.of.Spots.Under.Tissue")
# vars_number_pre <- vars_number[seq_len(split_point - 1)]
# vars_number_seq <- vars_number[-seq_len(split_point - 1)]
# vars_number_comb <- expand.grid(pre = vars_number_pre, seq = vars_number_seq)
# 
# ## Plot pre-sequencing vs post-sequencing
# pdf(
#     file.path(
#         dir_plots,
#         "spaceranger_scatterplot_pre-sequencing_vs_sequencing_metrics.pdf"
#     ),
#     useDingbats = FALSE,
#     width = length(unique(merged_metrics$study_name_sheet)) * 2,
#     height = length(unique(merged_metrics$study_name_sheet)) * 2
# )
# for (i in seq_len(nrow(vars_number_comb))) {
#     # for (i in 1:3) {
#     set.seed(20220914)
#     p <-
#         ggscatter(
#             merged_metrics_both,
#             x = as.character(vars_number_comb$pre[i]),
#             y = as.character(vars_number_comb$seq[i]),
#             color = "sample_status",
#             palette = sample_status_colors[grep("Seq", names(sample_status_colors))],
#             facet.by = "study_name_sheet_short",
#             add = "reg.line",
#             label = "slide_serial_capture_area",
#             repel = TRUE,
#             font.label = list(size = 10),
#             # add.params = list(size = 3),
#             legend = "top",
#             ggtheme = theme_pubr(base_size = 20),
#             cor.coef = TRUE,
#             # Add correlation coefficient. see ?stat_cor
#             cor.coeff.args = list(
#                 method = "pearson",
#                 label.x.npc = "center",
#                 label.y.npc = "bottom"
#             )
#         )
#     print(p)
# }
# dev.off()
# 
# ## Plot post-sequencing vs pre-sequencing
# vars_number_comb_rev <- expand.grid(seq = vars_number_seq, pre = vars_number_pre)
# pdf(
#     file.path(
#         dir_plots,
#         "spaceranger_scatterplot_sequencing_vs_pre-sequencing_metrics.pdf"
#     ),
#     useDingbats = FALSE,
#     width = length(unique(merged_metrics$study_name_sheet)) * 2,
#     height = length(unique(merged_metrics$study_name_sheet)) * 2
# )
# for (i in seq_len(nrow(vars_number_comb_rev))) {
#     # for (i in 1:3) {
#     set.seed(20220914)
#     p <-
#         ggscatter(
#             merged_metrics_both,
#             x = as.character(vars_number_comb_rev$pre[i]),
#             y = as.character(vars_number_comb_rev$seq[i]),
#             color = "sample_status",
#             palette = sample_status_colors[grep("Seq", names(sample_status_colors))],
#             facet.by = "study_name_sheet_short",
#             add = "reg.line",
#             label = "slide_serial_capture_area",
#             repel = TRUE,
#             font.label = list(size = 10),
#             # add.params = list(size = 3),
#             legend = "top",
#             ggtheme = theme_pubr(base_size = 20),
#             cor.coef = TRUE,
#             # Add correlation coefficient. see ?stat_cor
#             cor.coeff.args = list(
#                 method = "pearson",
#                 label.x.npc = "center",
#                 label.y.npc = "bottom"
#             )
#         )
#     print(p)
# }
# dev.off()
# 
# ## Plot pre-sequencing vs pre-sequencing
# vars_number_comb_pre <- data.frame(t(combn(vars_number_pre, 2)))
# colnames(vars_number_comb_pre) <- c("pre_var1", "pre_var2")
# pdf(
#     file.path(
#         dir_plots,
#         "spaceranger_scatterplot_pre-sequencing_vs_pre-sequencing_metrics.pdf"
#     ),
#     useDingbats = FALSE,
#     width = length(unique(merged_metrics$study_name_sheet)) * 2,
#     height = length(unique(merged_metrics$study_name_sheet)) * 2
# )
# for (i in seq_len(nrow(vars_number_comb_pre))) {
#     # for (i in 1:3) {
#     set.seed(20220914)
#     p <-
#         ggscatter(
#             subset(merged_metrics, !is.na(sheet_file)),
#             x = as.character(vars_number_comb_pre$pre_var1[i]),
#             y = as.character(vars_number_comb_pre$pre_var2[i]),
#             color = "sample_status",
#             palette = sample_status_colors,
#             facet.by = "study_name_sheet_short",
#             add = "reg.line",
#             label = "slide_serial_capture_area",
#             repel = TRUE,
#             font.label = list(size = 10),
#             # add.params = list(size = 3),
#             legend = "top",
#             ggtheme = theme_pubr(base_size = 20),
#             cor.coef = TRUE,
#             # Add correlation coefficient. see ?stat_cor
#             cor.coeff.args = list(
#                 method = "pearson",
#                 label.x.npc = "center",
#                 label.y.npc = "bottom"
#             )
#         )
#     print(p)
# }
# dev.off()
# 
# 
# ## Plot post-sequencing vs post-sequencing
# vars_number_comb_seq <- data.frame(t(combn(vars_number_seq, 2)))
# colnames(vars_number_comb_seq) <- c("seq_var1", "seq_var2")
# pdf(
#     file.path(
#         dir_plots,
#         "spaceranger_scatterplot_sequencing_vs_sequencing_metrics.pdf"
#     ),
#     useDingbats = FALSE,
#     width = length(unique(merged_metrics$study_name_sheet)) * 2,
#     height = length(unique(merged_metrics$study_name_sheet)) * 2
# )
# for (i in seq_len(nrow(vars_number_comb_seq))) {
#     # for (i in 1:3) {
#     set.seed(20220914)
#     p <-
#         ggscatter(
#             subset(merged_metrics, !is.na(summary_file)),
#             x = as.character(vars_number_comb_seq$seq_var1[i]),
#             y = as.character(vars_number_comb_seq$seq_var2[i]),
#             color = "sample_status",
#             palette = sample_status_colors[grep("Seq", names(sample_status_colors))],
#             facet.by = "study_name_sheet_short",
#             add = "reg.line",
#             label = "slide_serial_capture_area",
#             repel = TRUE,
#             font.label = list(size = 10),
#             # add.params = list(size = 3),
#             legend = "top",
#             ggtheme = theme_pubr(base_size = 20),
#             cor.coef = TRUE,
#             # Add correlation coefficient. see ?stat_cor
#             cor.coeff.args = list(
#                 method = "pearson",
#                 label.x.npc = "center",
#                 label.y.npc = "bottom"
#             )
#         )
#     print(p)
# }
# dev.off()
# 
# ## Reproducibility information
# print("Reproducibility information:")
# Sys.time()
# proc.time()
# options(width = 120)
# session_info()
