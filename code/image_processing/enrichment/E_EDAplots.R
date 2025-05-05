setwd("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/")

library("here")
library("ggplot2")
library(ggridges)
library(ggrepel)


Dr = here("processed-data", "image_processing")

data = read.csv(here(Dr, "AFseg", "data.csv"))

p = ggplot(data, x = Sample, y = Area) + geom_boxplot(aes(x = Sample, y = Area)) + ylim(0,50000) +
	geom_point(aes(x = Sample, y = Area), position = position_jitter(width = 0.1), size = 1, shape =1) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(here(Dr,"area.png"), plot = p, width = 8, height = 6, units = "in", dpi = 300)

png(here(Dr,"area.png"), width = 800, height = 600, units = "px", res = 300)
print(p)
dev.off()

p =ggplot(data, aes(x = Area, y = MeanIntensity, color = Sample)) +
   geom_point(alpha = 0.5, size = 1) + theme(legend.position = "none") +
   xlim(0, 50000) +
   geom_text_repel(
       data = subset(data, MeanIntensity < 1 & Area < 10),  # Subset data based on conditions
       aes(label = Sample),
       box.padding = 0.5,
       point.padding = 0.5
     ) 
   
p = ggplot(data, aes(x = Area, y = Sample, fill = Sample)) +
     geom_density_ridges(scale = 0.8, rel_min_height = 0.01) +
     theme_minimal()
	 
	 
filtered_data <- subset(data, MeanIntensity < 0.2 & Area < 500)
unique_samples <- unique(filtered_data$Sample)
unique_samples 



sample_counts <- table(data$Sample)

# Convert counts to data frame for plotting
sample_counts_df <- as.data.frame(sample_counts)
names(sample_counts_df) <- c("Sample", "Count")  # Rename columns

# Create bar plot
ggplot(sample_counts_df, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Occurrences of Each Sample",
       x = "Sample", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels
  
  
  sample_area_sum <- aggregate(Area ~ Sample, data = data, FUN = sum)

  # Create bar plot
  ggplot(sample_area_sum, aes(x = Sample, y = Area)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Total Area for Each Sample",
         x = "Sample", y = "Total Area") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels
  
  
  
	sample_subset <- subset(data, Sample == "V13F27-295_A1")
	max_area <- max(sample_subset$Area)
  
  
	ggplot(data, aes(x = Area)) +
	  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
	  labs(title = "Histogram of Area",
	       x = "Area", y = "Frequency") + xlim(0,50000) + ylim(0,500000)
	  theme_minimal()
	  


	  min_intensity <- aggregate(MeanIntensity ~ Sample, data = data, FUN = min)

	  # Create ggplot
	  ggplot(min_intensity, aes(x = Sample, y = MeanIntensity)) +
	    geom_point() +  # Add points
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



agg_data <- aggregate(cbind(MeanIntensity, Area) ~ Sample, data = data, FUN = function(x) c(min = min(x), sum = sum(x), me = mean(x), ma = max(x)))
ggplot(agg_data , aes(x = Sample, y = Area[,"ma"])) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Max Area for Each Sample",
       x = "Sample") + ylim(0,100000) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels


plot(agg_data$Area[,"ma"], 1:64)

filtered_data <- subset(data, Area < 500)
ggplot(agg_data, aes(x = MeanIntensity[,"ma"], y = Area[,"ma"], label = Sample)) +
  geom_point() +
  geom_text(vjust = -0.5) + 
  geom_hline(yintercept = 30000, linetype = "dashed", color = "red")+
  scale_y_log10()+scale_x_log10()