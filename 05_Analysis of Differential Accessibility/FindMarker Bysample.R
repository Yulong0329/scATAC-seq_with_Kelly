.libPaths("/home1/yulongqi/R/x86_64-pc-linux-gnu-library/4.4/Signac/R")
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

Bysample_seurat <- readRDS("/home1/yulongqi/Normalization_bysample&allsample/Normalization_seurat_Bysample.rds")

Bysample_seurat$seurat_clusters <- as.character(Bysample_seurat$seurat_clusters)
distinct_clusters <- c("10", "6", "9", "7", "15", "8", "11", "3", "16", "12")

#Create a metadata column named distinct_cluster, retaining only the cluster you specified and the rest as NA
Bysample_seurat$distinct_cluster <- ifelse(Bysample_seurat$seurat_clusters %in% distinct_clusters,
                             Bysample_seurat$seurat_clusters, NA)

Idents(Bysample_seurat) <- "distinct_cluster"

#Subsets the target cluster
Bysample_seurat_subset <- subset(Bysample_seurat, idents = distinct_clusters)

pdf("/home1/yulongqi/Feature_Matrix/UMAP_plot_Bysample.pdf", width = 8, height = 8)
DimPlot(Bysample_seurat_subset, label = TRUE) + NoLegend()


markers_list <- list()
#For each cluster, execute FindMarkers() separately: 
#compare with other clusters to find the accessible areas of differences
# FindMarkers: cluster vs the rest
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

saveRDS(markers_list, file = "/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_bysample/all_markers_by_cluster.rds")


