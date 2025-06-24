library(GenomicRanges)

# Read all-sample narrowPeak from batch
fpath <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_allsample/sbatch_peaks/SeuratProject_peaks.narrowPeak"
peak_df_all <- read.table(file = fpath,
                          col.names = c("chr", "start", "end", "name", "score", "strand",
                                        "fold_change", "neg_log10pvalue_summit",
                                        "neg_log10qvalue_summit", "relative_summit_position"),
                          stringsAsFactors = FALSE)

# Read by-sample GRanges object
load("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_bysample/Combined_and_reduced_Granges/combined_peaks.RData")

# Convert by-sample GRanges to data.frame
peak_df_by <- data.frame(
  seqnames = as.character(seqnames(combined.peaks)),
  start = start(combined.peaks),
  end = end(combined.peaks)
)

# Convert all-sample peak_df_all into same format (also rename chr to seqnames)
peak_df_all <- data.frame(
  seqnames = as.character(peak_df_all$chr),
  start = peak_df_all$start,
  end = peak_df_all$end
)

# Compute overlaps: for each all-sample peak (query), how many by-sample peaks (reference) overlap
peak_df_all$overlaps <- sapply(1:nrow(peak_df_all), function(i){
  start <- peak_df_all$start[i]
  end <- peak_df_all$end[i]
  chr <- peak_df_all$seqnames[i]
  sum(
    peak_df_by$seqnames == chr & (
      (peak_df_by$start <= start & peak_df_by$end >= start & peak_df_by$end <= end) |
        (peak_df_by$start >= start & peak_df_by$start <= end & peak_df_by$end >= end) |
        (peak_df_by$start >= start & peak_df_by$end <= end) |
        (peak_df_by$start <= start & peak_df_by$end >= end)
    )
  )
})
saveRDS(peak_df_all,file="/Users/yulongqiu/Desktop/Biostats/Master_Thesis/peak_df_all_with_overlaps.rds")

# Overlap from by-sample (query) to all-sample (reference)
peak_df_by$overlaps <- sapply(1:nrow(peak_df_by), function(i){
  start <- peak_df_by$start[i]
  end <- peak_df_by$end[i]
  chr <- peak_df_by$seqnames[i]
  sum(
    peak_df_all$seqnames == chr & (
      (peak_df_all$start <= start & peak_df_all$end >= start & peak_df_all$end <= end) |  # left partial
        (peak_df_all$start >= start & peak_df_all$start <= end & peak_df_all$end >= end) |  # right partial
        (peak_df_all$start >= start & peak_df_all$end <= end)                         |     # fully within
        (peak_df_all$start <= start & peak_df_all$end >= end)                               # fully contains
    )
  )
})
saveRDS(peak_df_by,file="/Users/yulongqiu/Desktop/Biostats/Master_Thesis/peak_df_by_with_overlaps.rds")

