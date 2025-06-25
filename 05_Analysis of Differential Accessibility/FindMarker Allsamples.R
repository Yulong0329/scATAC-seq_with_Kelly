.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

Allsample_seurat <- readRDS("/home1/yulongqi/Normalization_bysample&allsample/Normalization_seurat_0403.rds")

Allsample_seurat$seurat_clusters <- as.character(Allsample_seurat$seurat_clusters)
distinct_clusters <- c("3", "5", "6", "7", "8", "9", "10", "13", "15", "46", "49", "55")

Allsample_seurat$distinct_cluster <- ifelse(Allsample_seurat$seurat_clusters %in% distinct_clusters,
                                            Allsample_seurat$seurat_clusters, NA)

Idents(Allsample_seurat) <- "distinct_cluster"

#Subsets the target cluster
Allsample_seurat_subset <- subset(Allsample_seurat, idents = distinct_clusters)

table(Idents(Allsample_seurat_subset)) 
pdf("/home1/yulongqi/Feature_Matrix/UMAP_plot_allsample.pdf", width = 8, height = 8)
DimPlot(Allsample_seurat_subset, label = TRUE) + NoLegend()

markers_list <- list()
#For each cluster, execute FindMarkers() separately: 
#compare with other clusters to find the accessible areas of differences
for (cluster in distinct_clusters) {
  markers <- FindMarkers(
    a_subset,
    ident.1 = cluster,
    ident.2 = NULL,
    only.pos = FALSE,
    min.pct = 0,                     
    logfc.threshold = 0,            
    max.cells.per.ident = 5000      
  )
markers_list[[cluster]] <- markers
  saveRDS(markers, file = paste0("/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_bysample_0506/markers_cluster", cluster, "_vs_rest.rds"))
}

saveRDS(markers_list, file = "/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_allsample/all_markers_by_cluster.rds")


