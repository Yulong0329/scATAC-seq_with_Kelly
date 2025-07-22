.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(mclust)
Allsample_seurat <- readRDS("/home1/yulongqi/Normalization_bysample&allsample/Normalization_seurat_0403.rds")


Bysample_seurat <- readRDS("/home1/yulongqi/Normalization_bysample&allsample/Normalization_seurat_Bysample.rds")



Allsample_seurat$seurat_clusters <- as.character(Allsample_seurat$seurat_clusters)
Bysample_seurat$seurat_clusters <- as.character(Bysample_seurat$seurat_clusters)

# Find common cell
common_cells <- intersect(colnames(Allsample_seurat), colnames(Bysample_seurat))

# Build a dataframe and extract the cluster label from the two objects
cluster_df <- data.frame(
  cell = common_cells,
  allsample = Allsample_seurat$seurat_clusters[common_cells],
  bysample  = Bysample_seurat$seurat_clusters[common_cells]
)

# Calculation standard ARI (all shared cells, no "other")
ari_standard <- adjustedRandIndex(cluster_df$allsample, cluster_df$bysample)