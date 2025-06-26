library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ggplot2)
library(scales)  
library(rtracklayer)
library(patchwork)


######Promoter
# Load all-samples peak data
# Read all-sample narrowPeak from batch
fpath <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_allsample/sbatch_peaks/SeuratProject_peaks.narrowPeak"
peak_df_all <- read.table(file = fpath,
                          col.names = c("chr", "start", "end", "name", "score", "strand",
                                        "fold_change", "neg_log10pvalue_summit",
                                        "neg_log10qvalue_summit", "relative_summit_position"),
                          stringsAsFactors = FALSE)
# Convert each file to GRanges object
all_sample_granges <- makeGRangesFromDataFrame(df = peak_df_all, keep.extra.columns = TRUE)
length(all_sample_granges)

#operate by_Sample_Peak data
load("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_bysample/Combined_and_reduced_Granges/combined_peaks.RData")
combined.peaks


#Calculate Promoters in all sample and by sample
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters_data <- promoters(txdb, upstream = 2000, downstream = 200)

#Find overlap between promoters and by-sample peaks
overlaps_by_sample <- findOverlaps(promoters_data, combined.peaks)
overlaps_by_sample
#The number of unique promoters covered by the by-sample strategy
captured_promoters_reduced <- length(unique(queryHits(overlaps_by_sample)))
captured_promoters_reduced

#Find overlap between promoters and all-samples peaks
overlaps_all_sample <- findOverlaps(promoters_data, all_sample_granges)
#The number of unique promoters covered by the all-samples strategy
captured_promoters_all <- length(unique(queryHits(overlaps_all_sample)))
captured_promoters_all


###Check how many peaks overlap with promoters
# By Sample
peaks_with_promoters_by_sample <- length(unique(subjectHits(overlaps_by_sample)))
# All Sample
peaks_with_promoters_all_sample <- length(unique(subjectHits(overlaps_all_sample)))


# calculate percentages level
percent_by_sample <- round((peaks_with_promoters_by_sample / length(combined.peaks)) * 100, 1)
percent_all_sample <- round((peaks_with_promoters_all_sample / length(all_sample_granges)) * 100, 1)

#create df
df_peaks <- data.frame(
  Method = rep(c("By Sample", "All Samples"), each = 2),
  Type = rep(c("Overlap promoter", "Other peaks"), 2),
  Count = c(peaks_with_promoters_by_sample,
            length(combined.peaks) - peaks_with_promoters_by_sample,
            peaks_with_promoters_all_sample,
            length(all_sample_granges) - peaks_with_promoters_all_sample),
  Label = c(paste0(percent_by_sample, "%"), "", paste0(percent_all_sample, "%"), "")
)



# Define Colors
fill_colors <- c("Overlap promoter (By Sample)" = "salmon",
                 "Overlap promoter (All Samples)" = "skyblue",
                 "Other peaks" = "white")

# Set the "fill group" column Settings
df_peaks$fill_group <- with(df_peaks, ifelse(Method == "By Sample" & Type == "Overlap promoter",
                                             "Overlap promoter (By Sample)",
                                             ifelse(Method == "All Samples" & Type == "Overlap promoter",
                                                    "Overlap promoter (All Samples)",
                                                    "Other peaks")))

# Creating plot
plot_promoter<-ggplot(df_peaks, aes(x = Method, y = Count, fill = fill_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = fill_colors,
    name = NULL,
    labels = c(
      "Other peaks",
      "All-Samples peaks overlapping promotor",
      "By-Sample peaks overlapping promotor"
    )
  ) +
  scale_y_continuous(labels = comma) +
  labs(title = "What % of peaks overlap a promoter?",
       x = NULL,
       y = "Peak Count") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.position = "right")  
ggsave("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/promoter_peak_overlap.png",plot_promoter,
       width = 8, height = 6, dpi = 300)



#####Enhancer
enhancer_path <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Enhancer&Promoters/enhancer/ENCFF231VWX.bed.gz"

#read .bed
bed_df <- read.table(enhancer_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

head(bed_df)
tail(bed_df)

#convert into granges
enhancers <- GRanges(
       seqnames = bed_df[[1]],
       ranges = IRanges(start = bed_df[[2]] + 1, end = bed_df[[3]])  # BED 是 0-based
)
enhancers

#Calculate the number of enhancer that fall on all sample peaks
overlaps_all_sample_enhancer <- findOverlaps(enhancers, all_sample_granges)
captured_enhancers_all <- length(unique(queryHits(overlaps_all_sample_enhancer)))
cat("Captured enhancers (all sample peaks):", captured_enhancers_all, "\n")

#Calculate the number of enhancer that fall on by sample peaks
overlaps_reduced_enhancer <- findOverlaps(enhancers, combined.peaks)
captured_enhancers_reduced <- length(unique(queryHits(overlaps_reduced_enhancer)))
cat("Captured enhancers (by sample peaks):", captured_enhancers_reduced, "\n")




### Calculate how many peaks are overlapped with enhancer
# By Sample
peaks_with_enhancer_by_sample <- length(unique(subjectHits(overlaps_reduced_enhancer))); peaks_with_enhancer_by_sample
# All Sample
peaks_with_enhancer_all_sample <- length(unique(subjectHits(overlaps_all_sample_enhancer))); peaks_with_enhancer_all_sample


###Make plot
percent_by_sample <- round((peaks_with_enhancer_by_sample / length(combined.peaks)) * 100, 1)
percent_all_sample <- round((peaks_with_enhancer_all_sample / length(all_sample_granges)) * 100, 1)

df_enhancer <- data.frame(
  Method = rep(c("By Sample", "All Samples"), each = 2),
  Type = rep(c("Overlap enhancer", "Other peaks"), 2),
  Count = c(peaks_with_enhancer_by_sample,
            length(combined.peaks) - peaks_with_enhancer_by_sample,
            peaks_with_enhancer_all_sample,
            length(all_sample_granges) - peaks_with_enhancer_all_sample),
  Label = c(paste0(percent_by_sample, "%"), "", paste0(percent_all_sample, "%"), "")
)

###Creating Plot
fill_colors <- c("Other peaks" = "white",
                 "Overlap enhancer (All Samples)" = "skyblue",
                 "Overlap enhancer (By Sample)" = "salmon")

# Add the "fill_group" column for coloring
df_enhancer$fill_group <- with(df_enhancer, ifelse(Method == "By Sample" & Type == "Overlap enhancer",
                                                   "Overlap enhancer (By Sample)",
                                                   ifelse(Method == "All Samples" & Type == "Overlap enhancer",
                                                          "Overlap enhancer (All Samples)",
                                                          "Other peaks")))

# Creating Plot
plot_enhancer<-ggplot(df_enhancer, aes(x = Method, y = Count, fill = fill_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = fill_colors,
    name = NULL,
    labels = c(
      "Other peaks",
      "All-Samples peaks overlapping enhancer",
      "by-Sample peaks overlapping enhancer"
    )
  ) +
  scale_y_continuous(labels = comma) +  
  labs(title = "What % of peaks overlap an enhancer?",
       x = NULL,
       y = "Peak Count") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_blank())
ggsave("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/enhancer_peak_overlap.png",plot_enhancer,
       width = 8, height = 6, dpi = 300)



#######Combine Enhancers & Promoters and find peaks they all have
#Merge enhancers and promoters to study the overlap of the number of them combined with the two peaks
Enhancers_Promoters_ranges <- c(promoters_data, enhancers); Enhancers_Promoters_ranges
length(Enhancers_Promoters_ranges)

###Find overlaps between allsample peaks and bysample peaks
overlaps_Enhancers_Promoters_allsample <- findOverlaps(Enhancers_Promoters_ranges, all_sample_granges)
overlaps_Enhancers_Promoters_allsample

overlaps_Enhancers_Promoters_bysample <- findOverlaps(Enhancers_Promoters_ranges, combined.peaks)
overlaps_Enhancers_Promoters_bysample


###Peaks overlaps with either Enhancers or Promoters
peaks_with_Enhancers_Promoters_by_sample <- length(unique(subjectHits(overlaps_Enhancers_Promoters_bysample)))
# All Sample
peaks_with_Enhancers_Promoters_all_sample <- length(unique(subjectHits(overlaps_Enhancers_Promoters_allsample)))
# print
peaks_with_Enhancers_Promoters_by_sample
peaks_with_Enhancers_Promoters_all_sample


###Make plot
percent_by_sample <- round((peaks_with_Enhancers_Promoters_by_sample / length(combined.peaks)) * 100, 1)
percent_all_sample <- round((peaks_with_Enhancers_Promoters_all_sample / length(all_sample_granges)) * 100, 1)

df_enhancer <- data.frame(
  Method = rep(c("By Sample", "All Samples"), each = 2),
  Type = rep(c("Overlap enhancer", "Other peaks"), 2),
  Count = c(peaks_with_Enhancers_Promoters_by_sample,
            length(combined.peaks) - peaks_with_Enhancers_Promoters_by_sample,
            peaks_with_Enhancers_Promoters_all_sample,
            length(all_sample_granges) - peaks_with_Enhancers_Promoters_all_sample),
  Label = c(paste0(percent_by_sample, "%"), "", paste0(percent_all_sample, "%"), "")
)

###Recreate the merged enhancer and promoter
fill_colors <- c("Other peaks" = "white",
                 "Overlap enhancer (All Samples)" = "skyblue",
                 "Overlap enhancer (By Sample)" = "salmon")
df_enhancer$fill_group <- with(df_enhancer, ifelse(Method == "By Sample" & Type == "Overlap enhancer",
                                                   "Overlap enhancer (By Sample)",
                                                   ifelse(Method == "All Samples" & Type == "Overlap enhancer",
                                                          "Overlap enhancer (All Samples)",
                                                          "Other peaks")))

# Creating Plot
plot_both<-ggplot(df_enhancer, aes(x = Method, y = Count, fill = fill_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = fill_colors,
    name = NULL,
    labels = c(
      "Other peaks",
      "All-samples peaks overlapping either promoter or enhancer",
      "By-Sample peaks overlapping either promoter or enhancer"
    )
  ) +
  scale_y_continuous(labels = comma) +  
  labs(title = "What % of peaks overlap either a enhancer or a promoter?",
       x = NULL,
       y = "Peak Count") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "right")
ggsave("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/Enhancers_Promoters_peak_overlap.png",plot_both,
       width = 8, height = 6, dpi = 300)

###Mark ABC and display all the graphs as two in the first row and one in the second row
combined_plot <- (
  (plot_promoter | plot_enhancer) / 
    (plot_both     | plot_spacer())
) +
  plot_layout(widths = c(1, 0.4)) +
  plot_annotation(tag_levels = 'A')  # Automatically add A, B, and C







