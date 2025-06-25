library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

Bysample_peaks <- readRDS('/home1/yulongqi/Feature_Matrix/FeatureMatrix_Bysample.rds')
Bysample_peaks <- as(Bysample_peaks, "dgCMatrix")

###Create a Seurat object
# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(counts = Bysample_peaks)
# Create Seurat object
Bysample_seurat <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

### Normalization and linear dimensional reduction
Bysample_seurat <- RunTFIDF(Bysample_seurat)
Bysample_seurat <- FindTopFeatures(Bysample_seurat, min.cutoff = 'q0')
Bysample_seurat <- RunSVD(Bysample_seurat)

### Assess the correlation between each LSI component and sequencing depth
depth_plot <- DepthCor(Bysample_seurat)
ggsave(filename = "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/LSI_vs_Depth_correlation.png", plot = depth_plot, width = 6, height = 4, dpi = 300)

### Non-linear dimension reduction and clustering
Bysample_seurat <- RunUMAP(object = Bysample_seurat, reduction = 'lsi', dims = 2:30)
Bysample_seurat <- FindNeighbors(object = Bysample_seurat, reduction = 'lsi', dims = 2:30)
Bysample_seurat <- FindClusters(object = Bysample_seurat, verbose = FALSE, algorithm = 3)

### Save dimensionally reduced and clustered Seurat object
saveRDS(Bysample_seurat, file = "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Bysample_seurat.rds")

### UMAP plot and save
p <- DimPlot(object = Bysample_seurat, label = TRUE) + NoLegend()
ggsave(filename = "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/UMAP_Bysample_clusters.png", plot = p, width = 6, height = 5, dpi = 300)






