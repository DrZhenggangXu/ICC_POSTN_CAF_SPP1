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

sce_malignant <- subset(sce,cell_type == 'Malignant')
sce_malignant <- NormalizeData(sce_malignant,
                               normalization.method = 'LogNormalize',
                               scale.factor = 1e4)
sce_malignant <- FindVariableFeatures(sce_malignant,
                                      selection.method = "vst",
                                      nfeatures = 2000)
sce_malignant <- ScaleData(sce_malignant,
                           vars.to.regress = c("nCount_RNA","percent.mt"))
sce_malignant <- RunPCA(sce_malignant,features = VariableFeatures(sce_malignant))
sce_malignant <- RunHarmony(sce_malignant,
                      group.by.vars = 'orig.ident'
                      )
dims <- 10
sce_malignant <- RunUMAP(sce_malignant, reduction = "harmony", dims = 1:dims)
sce_malignant <- RunTSNE(sce_malignant, reduction = "harmony", dims = 1:dims)
sce_malignant <- FindNeighbors(sce_malignant, reduction = "harmony", dims = 1:dims)
res <- 0.1
sce_malignant <- FindClusters(sce_malignant,resolution = res)
sce_malignant$cell_type <- paste0('Tumor_c',sce_malignant$seurat_clusters)
DimPlot(sce_malignant,reduction="umap",
        group.by="cell_type",
        shuffle = F,
        pt.size=1.0,raster=FALSE,
        label=F,label.box=F,label.size=6)
FeaturePlot(sce_malignant,features = 'SPP1',reduction = 'umap',raster = F,
            pt.size = 1.0,
            cols = c('grey','red'))
sce_malignant@meta.data[sce_malignant$cell_type=="Tumor_c0",
                        'cell_type'] <- 'SPP1+Malignant'
sce_malignant@meta.data[sce_malignant$cell_type=="Tumor_c1",
                        'cell_type'] <- 'SPP1-Malignant'
sce_malignant@meta.data[sce_malignant$cell_type=="Tumor_c2",
                        'cell_type'] <- 'SPP1+Malignant'
sce_malignant@meta.data[sce_malignant$cell_type=="Tumor_c3",
                        'cell_type'] <- 'SPP1-Malignant'
sce_malignant@meta.data[sce_malignant$cell_type=="Tumor_c4",
                        'cell_type'] <- 'SPP1+Malignant'
sce_malignant@meta.data[sce_malignant$cell_type=="Tumor_c5",
                        'cell_type'] <- 'SPP1+Malignant'
DimPlot(sce_malignant,reduction="umap",
        group.by="cell_type",
        pt.size=1.0)
Idents(sce_malignant) <- 'cell_type'
markers <- FindAllMarkers(sce_malignant,
                          logfc.threshold = 0.25,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.25,
                          only.pos = F
                          )
jjVolcano(diffData = markers,size = 4.0,
          topGeneN = 5,
          fontface = 'italic')

