getwd()
library(Seurat)
library(limma)
library(dplyr)
library(tidyverse)
library(patchwork)
library(devtools)
library(harmony)
library(DoubletFinder)
library(multtest)
library(Rcpp)
library(ggplot2)
library(clustree)
library(org.Hs.eg.db)
library(data.table)
library(rtracklayer)
library(R.utils)
library(ggpubr)

sce_fibro <- subset(sce,cell_type == 'Fibroblasts')
sce_fibro <- NormalizeData(sce_fibro,
                           normalization.method = 'LogNormalize',
                           scale.factor = 1e4)
sce_fibro <- FindVariableFeatures(sce_fibro,
                                  selection.method = "vst",
                                  nfeatures = 2000)
sce_fibro <- ScaleData(sce_fibro,
                       vars.to.regress = c("nCount_RNA","percent.mt"))
sce_fibro <- RunPCA(sce_fibro,features = VariableFeatures(sce_fibro))
sce_fibro <- RunHarmony(sce_fibro,
                        group.by.vars = 'orig.ident')
dims <- 10
sce_fibro <- RunUMAP(sce_fibro, reduction = "harmony", dims = 1:dims)
sce_fibro <- RunTSNE(sce_fibro, reduction = "harmony", dims = 1:dims)
sce_fibro <- FindNeighbors(sce_fibro, reduction = "harmony", dims = 1:dims)
res <- 0.2
sce_fibro <- FindClusters(sce_fibro,resolution = res)
markers <- FindAllMarkers(sce_fibro,
                          logfc.threshold = 0.25,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.1,
                          only.pos = TRUE)
sce_fibro@meta.data[sce_fibro$seurat_clusters==0,
                    'cell_type'] <- 'CAF_c0_RGS5'
sce_fibro@meta.data[sce_fibro$seurat_clusters==1,
                    'cell_type'] <- 'CAF_c1_POSTN'
sce_fibro@meta.data[sce_fibro$seurat_clusters==2,
                    'cell_type'] <- 'CAF_c2_CD74'
sce_fibro@meta.data[sce_fibro$seurat_clusters==3,
                    'cell_type'] <- 'CAF_c3_APOD'
DimPlot(sce_fibro,reduction="umap",
        group.by="cell_type")
markers_1 <- FindAllMarkers(sce_fibro,
                          logfc.threshold = 0.5,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.4,
                          only.pos = F)
jjVolcano(diffData = markers_1,
          log2FC.cutoff = 0.4,
          size = 4.0,
          topGeneN = 5,
          fontface = 'italic')+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 15, colour = 'black', angle = 45, hjust = 1, vjust = .5), 
    legend.text = element_text(size = 14),
    legend.position = "top"
  )
