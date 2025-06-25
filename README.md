# scATAC-seq_with_Kelly
# About
This repository contains the full analysis code and data processing scripts for my master's thesis project, which compares two peak calling strategies—by-sample and all-sample—in single-cell ATAC-seq data from colorectal cancer samples.

## Project Overview
The project systematically compares the impact of two peak calling strategies on:

Peak length distribution
Overlap with known regulatory elements (promoters and enhancers)
Peaks Overlap Comparison
Downstream analysis outcomes (Normalization, Dimensional Reduction, Clustering, Differential Accessibility)

# Tools & Packages
R (v4.4.1)
Key packages: Signac, Seurat, ggplot2, GenomicRanges, TxDb.Hsapiens.UCSC.hg38.knownGene, ENCODE(H3K27ac)
Peak calling: MACS3, run via Signac or batch script

# Notes
Peak calling was performed on 89 colorectal cancer samples
Only chromosome 1 was used in downstream analyses due to memory limitations
