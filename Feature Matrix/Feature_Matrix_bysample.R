.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)

fpaths <- system("ls /Users/yulongqiu/Desktop/Biostats/Master_Thesis/Sample/*.tsv.gz", intern = TRUE)

fragments <- lapply(fpaths, function(fp){
    total_counts <- CountFragments(fp)#Count how many fragments appear in each barcode (i.e., a single cell).
    barcodes <- total_counts$CB[total_counts$frequency_count > 100]
    CreateFragmentObject(fp, cells = barcodes)#Store the meta-information related to the fragments
})

#Create bysample peaks subset(chr1)
features <- readRDS('/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_bysample/Combined_and_reduced_Granges/combined_peaks.rds')

#Create subset of by sample peaks(only chr1 in the dataset)
subset_chr1_bysample <- subset(features, seqnames == "chr1")
subset_chr1_bysample

# Generate the feature matrix
feature_matrix <- FeatureMatrix(
  fragments = fragments,
  features = subset_chr1_bysample
)

saveRDS(feature_matrix, file = "/home1/yulongqi/Feature_Matrix/FeatureMatrix_bysample.rds")




