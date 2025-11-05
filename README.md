# scRNA-seq-Analysis-Seurat-BF528
This repository stores the code for a scRNA-seq pipeline built to complete a final project for Boston University's course BF528. The pipeline includes quality control, pre-processing, and seurat workflow for single cell analysis. Additionally, the workflow featues single cell cluster annotation using both automated and manual methods. 

The data used for this project was pulled from the publiation 'Identification of distinct tumor cell populations and key genetic mechanisms through single cell sequencing in hepatoblastoma' produced by Bondoc et. al. (Bondoc et. al., 2021). 

Two figures were reproduced from the intial paper in order to validate the pieplines results. Figure 3b was reproduced from the paper and stored in this repository as 'marker_umap_grid.png'. It is important to note that the figure is not an exact carbon copy from the original paper, but this is expected as the workflow analysis are different. Additionally, the figure produced by this workflow matches expression patterns of the biomarkers in respective cell types. 

## Requirements

R (>= 4.0 recommended)

Packages: Seurat, SeuratWrappers, tidyverse, ggplot2, DoubletFinder, devtools, SeuratWrappers, Harmony (or HarmonyIntegration wrapper via SeuratWrappers), presto (installed from GitHub), SingleR, celldex, SummarizedExperiment, patchwork

## Inputs

Seven 10X filtered_feature_bc_matrix or raw_feature_bc_matrix directories (one per sample) named as in the notebook (e.g. HB17_background_raw_feature_bc_matrix).

Workflow (high-level)

Import — Read10X and CreateSeuratObject for each sample (min.cells = 3, min.features = 200).

QC — compute percent.mt (MT-), inspect violin plots for nFeature_RNA, nCount_RNA, percent.mt and set per-sample thresholds.

Filtering — subset cells by nFeature_RNA, nCount_RNA, and percent.mt (per-sample thresholds chosen from violin plots).

Per-sample analysis — NormalizeData, FindVariableFeatures, ScaleData, RunPCA (used 20 PCs), FindNeighbors, FindClusters (resolution = 0.5), RunUMAP.

Doublet detection — paramSweep / summarizeSweep → find.pK → estimate expected doublets → doubletFinder (prefix DD_ added).

Merge — merge DoubletFinder-corrected Seurat objects into a single Seurat object.

Integrated analysis — Normalize, FindVariableFeatures, ScaleData, RunPCA, ElbowPlot → FindNeighbors/FindClusters/RunUMAP.

Integration — IntegrateLayers with Harmony (creates harmony reduction); re-cluster and run UMAP on harmony embedding.

Marker discovery — JoinLayers, set identities to harmony clusters, FindAllMarkers (only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), extract top 5 markers per cluster.

Annotation — Automated annotation with SingleR (celldex HumanPrimaryCellAtlasData) + manual curation from marker genes (cellxgene reference comparisons).

Reproduced figures — FeaturePlot grid (saved as marker_umap_grid.png), correlation scatterplots comparing mean-normalized expression between sample groups (report R² from linear models).

## Outputs

marker_umap_grid.png — combined FeaturePlots for selected genes

R objects in RMarkdown environment (Seurat objects at multiple stages)

Tables: markers, top5_markers

UMAP and correlation figures reproduced for the paper

Work Cited

1. Bondoc, A., Glaser, K., Jin, K. et al. Identification of distinct tumor cell populations and key genetic mechanisms through single cell
sequencing in hepatoblastoma. Commun Biol 4, 1049 (2021). https://doi.org/10.1038/s42003-021-02562-8
