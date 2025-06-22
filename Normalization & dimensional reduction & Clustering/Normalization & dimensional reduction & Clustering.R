library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

Allsample_peaks <- readRDS('/home1/yulongqi/Feature_Matrix/FeatureMatrix_allsample_0402.rds')
Allsample_peaks <- as(Allsample_peaks, "dgCMatrix")

###Create a Seurat object
# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(counts = Allsample_peaks)
# Create Seurat object
allsample_seurat <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

### Normalization and linear dimensional reduction
allsample_seurat <- RunTFIDF(allsample_seurat)
allsample_seurat <- FindTopFeatures(allsample_seurat, min.cutoff = 'q0')
allsample_seurat <- RunSVD(allsample_seurat)

### Assess the correlation between each LSI component and sequencing depth
depth_plot <- DepthCor(allsample_seurat)
ggsave(filename = "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/LSI_vs_Depth_correlation.png", plot = depth_plot, width = 6, height = 4, dpi = 300)

### Non-linear dimension reduction and clustering
allsample_seurat <- RunUMAP(object = allsample_seurat, reduction = 'lsi', dims = 2:30)
allsample_seurat <- FindNeighbors(object = allsample_seurat, reduction = 'lsi', dims = 2:30)
allsample_seurat <- FindClusters(object = allsample_seurat, verbose = FALSE, algorithm = 3)

### Save dimensionally reduced and clustered Seurat object
saveRDS(allsample_seurat, file = "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Allsample_seurat.rds")

### UMAP plot and save
p <- DimPlot(object = allsample_seurat, label = TRUE) + NoLegend()
ggsave(filename = "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/UMAP_allsample_clusters.png", plot = p, width = 6, height = 5, dpi = 300)






