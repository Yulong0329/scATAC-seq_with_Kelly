library(GenomicRanges)
library(dplyr)
library(ggplot2)

#load by_sample peak data
load("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_bysample/Combined_and_reduced_Granges/combined_peaks.RData")

# Load all-samples peak data
# Read all-sample narrowPeak from batch
fpath <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_allsample/sbatch_peaks/SeuratProject_peaks.narrowPeak"
peak_df_all <- read.table(file = fpath,
                          col.names = c("chr", "start", "end", "name", "score", "strand",
                                        "fold_change", "neg_log10pvalue_summit",
                                        "neg_log10qvalue_summit", "relative_summit_position"),
                          stringsAsFactors = FALSE)

# Convert all-sample peak_df_all into same format (also rename chr to seqnames)
all_sample_df <- data.frame(
  seqnames = as.character(peak_df_all$chr),
  start = peak_df_all$start,
  end = peak_df_all$end
)

# Convert each file to GRanges object
all_sample_granges <- makeGRangesFromDataFrame(df = all_sample_df, keep.extra.columns = TRUE)


# Create data frame with length and strategy labels
all_sample_df <- data.frame(length = width(all_sample_granges), strategy = "All Sample")
by_sample_df  <- data.frame(length = width(combined.peaks), strategy = "By Sample")
peak_data <- rbind(all_sample_df, by_sample_df)


# Set bin width and assign bins
bin_width <- 100
peak_data <- peak_data %>%
  mutate(bin_start = floor(length / bin_width) * bin_width)#Allocate each length to each bin

# Compute percentage for each bin
peak_summary <- peak_data %>%
  group_by(strategy, bin_start) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(strategy) %>%
  mutate(percent = count / sum(count) * 100)


###Creating plot of Length Distribution
#Define the output folder for saving plots
output_folder <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output"

peak_summary$strategy <- factor(peak_summary$strategy, levels = c("All Sample", "By Sample"))

p_combined <- ggplot(peak_summary, aes(x = bin_start, y = percent, fill = strategy)) +
  geom_col(alpha = 0.8, width = bin_width * 0.9, show.legend = TRUE) +
  facet_grid(. ~ strategy) +  
  scale_fill_manual(
    name = "Strategy",
    values = c("All Sample" = "skyblue", "By Sample" = "salmon"),
    labels = c("All-samples method", "By-sample method")
  ) +
  scale_x_continuous(
    name = "Peak Length (bp)",
    breaks = seq(0, 5000, by = 1000)
  ) +
  ylab("Percentage of Peaks (%)") +
  ggtitle("Peak Length Distribution by Strategy (Percentage)") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),  
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave(filename = file.path(output_folder, "peak_length_distribution_Percentage_combined2.png"), 
       plot = p_combined, width = 8, height = 6, dpi = 300)





