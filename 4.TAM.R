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
library(scRNAtoolVis)
options(stringsAsFactors = F)

sce_tam <- subset(sce,cell_type == 'TAMs')
sce_tam <- NormalizeData(sce_tam,
                         normalization.method = 'LogNormalize',
                         scale.factor = 1e4)
sce_tam <- FindVariableFeatures(sce_tam,
                                selection.method = "vst",
                                nfeatures = 2000)
sce_tam <- ScaleData(sce_tam,
                     vars.to.regress = c("nCount_RNA","percent.mt"))
sce_tam <- RunPCA(sce_tam,features = VariableFeatures(sce_tam))
sce_tam <- RunHarmony(sce_tam,
                      group.by.vars = 'orig.ident')
dims <- 10
sce_tam <- RunUMAP(sce_tam, reduction = "harmony", dims = 1:dims)
sce_tam <- RunTSNE(sce_tam, reduction = "harmony", dims = 1:dims)
sce_tam <- FindNeighbors(sce_tam, reduction = "harmony", dims = 1:dims)
res <- 0.1
sce_tam <- FindClusters(sce_tam,resolution = res)
markers <- FindAllMarkers(sce_tam,
                          logfc.threshold = 0.25,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.1,
                          only.pos = TRUE)
cluster2celltype <- c("0"="TAM_c0_SPP1",
                      "1"="TAM_c1_FCGR3B",
                      "2"="TAM_c2_THBS1", 
                      "3"="TAM_c3_IL8",
                      "4"="TAM_C4_HSPH1")
sce_tam[['cell_type']] = unname(cluster2celltype[sce_tam@meta.data$seurat_clusters])
DimPlot(sce_tam, reduction = 'umap', group.by = 'cell_type',
        label = F)
FeaturePlot(sce_tam,features = 'SPP1',reduction = 'umap',raster = F,
            cols = c('grey','red'))
Idents(sce_tam) <- 'cell_type'
markers <- FindAllMarkers(sce_tam,
                          logfc.threshold = 0.5,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.2,
                          only.pos = F)
jjVolcano(diffData = markers,
          size = 5.0,
          topGeneN = 5,
          fontface = 'italic')
