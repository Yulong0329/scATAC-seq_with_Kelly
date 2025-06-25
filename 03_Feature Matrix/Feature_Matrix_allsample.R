.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)

fpaths <- system("ls /project/kellystr_1320/GEO_Data/*fragments.tsv.gz", intern = TRUE)
fragments <- lapply(fpaths, function(fp){
    total_counts <- CountFragments(fp)
    barcodes <- total_counts$CB[total_counts$frequency_count > 100]
    CreateFragmentObject(fp, cells = barcodes)
})

# Define features
fpath <- "/project/kellystr_1320/GEO_Data/peaks_allsamples/SeuratProject_peaks.narrowPeak"

#features <- granges(fpath)
require(GenomicRanges)
df <- read.table(file = fpath,
                 col.names = c("chr", "start", "end", "name", "score", "strand",
                               "fold_change", "neg_log10pvalue_summit",
                               "neg_log10qvalue_summit", "relative_summit_position"))
# Convert each file to GRanges object and assign to 'features'
features <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)


#Create allsample peaks subset(chr1)
#Create subset of all sample peaks(only chr1 in the dataset)
features <- subset(features, seqnames == "chr1")


# Generate the feature matrix
feature_matrix <- FeatureMatrix(
  fragments = fragments,
  features = features
)

saveRDS(feature_matrix, file = "/home1/yulongqi/Feature_Matrix/FeatureMatrix_allsample_0402.rds")




