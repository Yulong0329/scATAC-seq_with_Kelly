.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)

fpaths <- system("ls /Users/yulongqiu/Desktop/Biostats/Master_Thesis/Sample/*.tsv.gz", intern = TRUE)
#intern = TRUE，让结果作为字符向量返回，而不是打印到终端
#system函数这一段列出所有的tsv.gz, 并把这些文件路径作为字符向量传回 R
fragments <- lapply(fpaths, function(fp){
    total_counts <- CountFragments(fp)
#这个函数是 Signac 包中的工具函数，用于读取 .fragments.tsv.gz 文件，
#统计每个 barcode（即单细胞）中出现了多少个 fragment（片段）
    barcodes <- total_counts$CB[total_counts$frequency_count > 100]
#这一步是筛选出片段数超过 100 的细胞，认为这些细胞可能是质量比较高的细胞。
    CreateFragmentObject(fp, cells = barcodes)
#这一步是创建一个 Fragment 对象—— Signac 中的一个数据结构，专门用来存储与片段有关的元信息。
#cells = barcodes 只保留那些你认为是“高质量”的细胞。
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




