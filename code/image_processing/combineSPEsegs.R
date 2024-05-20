
spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_wo_spg_N63.rds"))
PRECAST_df <- readRDS(here("processed-data/rds/spatial_cluster","PRECAST", "test_clus_label_df_semi_inform_k_2-16.rds"))
stopifnot(nrow(PRECAST_df) == ncol(spe))

precast_vars <- grep("^PRECAST_", colnames(PRECAST_df),value = TRUE)
spe <- spe[, spe$key %in% PRECAST_df$key]
# raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
col_data_df <- PRECAST_df |> right_join(
    colData(spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )

rownames(col_data_df) <- colnames(spe)
colData(spe) <- DataFrame(col_data_df)


library("readxl")
library("here")
library("tidyr")
library("spatialLIBD")
samples = read_excel(here("raw-data","experiment_info","VisiumSPG_PNN_Master.xlsx"))
samples <- as.data.frame(samples %>% select('Slide #', 'Array #'))
samples$sample_id = paste0(samples[,1],"_",samples[,2])
samples$SPpath = paste0('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/',samples$sample_id,'/outs/spatial/tissue_spot_counts.csv')
segmentations_list <-
  lapply(samples$sample_id, function(sample_id) {
    file <-samples$SPpath[samples$sample_id == sample_id]
    if (!file.exists(file)) {
      return(NULL)
    }
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sample_id)
    return(x)
  })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
  Reduce(function(...) {
    merge(..., all = TRUE)
  }, segmentations_list[lengths(segmentations_list) > 0])


## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spe) <- cbind(colData(spe), segmentation_info)


library("ggplot2")
coldata_df = as.data.frame(colData(spe))
coldata_df$IWFA <- ifelse(is.na(coldata_df$IWFA), 0, coldata_df$IWFA)
filtered_data <- coldata_df[coldata_df$IWFA < 1, ]


ggplot(coldata_df, aes(x = dx, y = PWFA)) + geom_boxplot() +
  labs(title = "PWFA by Dx",x = "Sex",y = "PWFA") + geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal()

ggplot(coldata_df, aes(x = sex, y = PWFA)) + geom_boxplot() +
    labs(title = "PWFA by Sex",x = "Sex",y = "PWFA") + geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal()
	
ggplot(coldata_df, aes(x = PWFA)) + geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
    labs(title = "Distribution of PWFA", x = "PWFA", y = "Frequency") +
    theme_minimal()

coldata_df_dx = coldata_df[coldata_df$dx == 'ntc', ]
proportion_data <- coldata_df_dx %>%
  group_by(sample_id, PRECAST_10, CNWFA_GT_0 = CNWFA > 0) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(sample_id, PRECAST_10) %>%
  mutate(Proportion = Count / sum(Count))

	# Show the first few rows of proportion_data
	head(proportion_data)

	# Create stacked bar plot
ggplot(proportion_data, aes(x = PRECAST_10, y = Proportion, fill = CNWFA_GT_0)) +
	  geom_bar(stat = "identity", position = "stack") +
	  labs(title = "Proportion of CNWFA>0 by Sample ID", x = "Sample ID", y = "Proportion", fill = "CNWFA > 0") +
	  theme_minimal()+facet_wrap(~ sample_id)
	  
coldata_df_dx = coldata_df[coldata_df$dx == 'scz', ]
proportion_data <- coldata_df_dx %>%
  group_by(sample_id, PRECAST_10, CNWFA_GT_0 = CNWFA > 0) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(sample_id, PRECAST_10) %>%
  mutate(Proportion = Count / sum(Count))

	  	# Show the first few rows of proportion_data
head(proportion_data)

	  	# Create stacked bar plot
ggplot(proportion_data, aes(x = PRECAST_10, y = Proportion, fill = CNWFA_GT_0)) +
	  geom_bar(stat = "identity", position = "stack") +
	  labs(title = "Proportion of CNWFA>0 by Sample ID", x = "Sample ID", y = "Proportion", fill = "CNWFA > 0") +
	  theme_minimal()+facet_wrap(~ sample_id)
	  
	  
	  coldata_df_dx = coldata_df[coldata_df$dx == 'ntc', ]
	  proportion_data <- coldata_df_dx %>%
	    group_by(sample_id, PRECAST_10, CNWFA_GT_0 = PWFA > 0.05) %>%
	    summarise(Count = n()) %>%
	    ungroup() %>%
	    group_by(sample_id, PRECAST_10) %>%
	    mutate(Proportion = Count / sum(Count))

	  	# Show the first few rows of proportion_data
	  	head(proportion_data)

	  	# Create stacked bar plot
	  ggplot(proportion_data, aes(x = PRECAST_10, y = Proportion, fill = CNWFA_GT_0)) +
	  	  geom_bar(stat = "identity", position = "stack") +
	  	  labs(title = "Proportion of CNWFA>0 by Sample ID", x = "Sample ID", y = "Proportion", fill = "CNWFA > 0") +
	  	  theme_minimal()+facet_wrap(~ sample_id)
	  
	  coldata_df_dx = coldata_df[coldata_df$dx == 'scz', ]
	  proportion_data <- coldata_df_dx %>%
	    group_by(sample_id, PRECAST_10, CNWFA_GT_0 = PWFA > 0.05) %>%
	    summarise(Count = n()) %>%
	    ungroup() %>%
	    group_by(sample_id, PRECAST_10) %>%
	    mutate(Proportion = Count / sum(Count))

	  	  	# Show the first few rows of proportion_data
	  head(proportion_data)

	  	  	# Create stacked bar plot
	  ggplot(proportion_data, aes(x = PRECAST_10, y = Proportion, fill = CNWFA_GT_0)) +
	  	  geom_bar(stat = "identity", position = "stack") +
	  	  labs(title = "Proportion of CNWFA>0 by Sample ID", x = "Sample ID", y = "Proportion", fill = "CNWFA > 0") +
	  	  theme_minimal()+facet_wrap(~ sample_id)