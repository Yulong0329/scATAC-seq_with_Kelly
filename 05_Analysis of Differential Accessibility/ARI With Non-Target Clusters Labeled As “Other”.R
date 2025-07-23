.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

Allsample_seurat <- readRDS("/home1/yulongqi/Normalization_bysample&allsample/Normalization_seurat_0403.rds")
Bysample_seurat <- readRDS("/home1/yulongqi/Normalization_bysample&allsample/Normalization_seurat_Bysample.rds")


Allsample_seurat$seurat_clusters <- as.character(Allsample_seurat$seurat_clusters)
Bysample_seurat$seurat_clusters <- as.character(Bysample_seurat$seurat_clusters)

distinct_clusters_all <- c("3", "5", "6", "7", "8", "9", "10", "13", "15", "46", "49", "55")
distinct_clusters_by <- c("10", "6", "9", "7", "15", "8", "11", "3", "16", "12")

Allsample_seurat$distinct_cluster <- ifelse(
  Allsample_seurat$seurat_clusters %in% distinct_clusters_all,
  Allsample_seurat$seurat_clusters,
  "other"
)

Bysample_seurat$distinct_cluster <- ifelse(
  Bysample_seurat$seurat_clusters %in% distinct_clusters_by,
  Bysample_seurat$seurat_clusters,
  "other"
)


common_cells <- intersect(colnames(Allsample_seurat), colnames(Bysample_seurat))

cluster_df <- data.frame(
  cell     = common_cells,
  allsample = Allsample_seurat[, common_cells]$distinct_cluster,
  bysample  = Bysample_seurat[, common_cells]$distinct_cluster)

adjustedRandIndex(cluster_df$allsample, cluster_df$bysample)