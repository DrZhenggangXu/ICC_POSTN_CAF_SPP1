getwd() 
rm(list = ls())
library(Seurat) # Version 4
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
library(ggsci)
options(stringsAsFactors = F)

## Data_Integration
sceList <- readRDS('raw_sceList.rds')
for (i in 1:length(sceList)){
  sce <- sceList[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sceList[[i]] <- sce
  }
sceList_filtered <- list()
for (i in 1:length(sceList)){
  sce <- sceList[[i]]
  sce <- subset(sce,
                nCount_RNA > 500 &
                  nCount_RNA < 8000 &
                  nFeature_RNA < 5000 &
                  nFeature_RNA > 300 &
                  percent.mt < 15 &
                  percent.hb < 5)
  sceList_filtered[[i]] <- sce
  }
sce <- merge(sceList_filtered[[1]],
             sceList_filtered[2:length(sceList_filtered)])
sce <- NormalizeData(object = sce,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4)
sce <- FindVariableFeatures(object = sce, 
                            selection.method = "vst", 
                            nfeatures = 2000)
sce <- ScaleData(object = sce,
                 vars.to.regress = c("nCount_RNA","percent.mt"))
sce <- RunPCA(sce)
sce <- RunHarmony(sce,group.by.vars = 'orig.ident')
dims <- 20
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:dims)
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:dims)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:dims)
res <- 0.3
sce <- FindClusters(sce,resolution = res) 
DimPlot(sce,reduction = 'umap',
        group.by = 'seurat_clusters')

## Annotation
genes_to_check <- c("EPCAM","KRT19","KRT7", # Malignant
                    "FYXD2", "TM4SF4", "ANXA4", #Cholangiocyte
                    "APOC3", "FABP1", "APOA1", # Hepatocye
                    "PTPRC", # immune
                    "CD19","CD79A","MS4A1","MZB1", #B
                    "CD2","CD3D","CD3E", #T
                    "CD4","IL2RA",
                    "CD8A","GZMB",'NKG7',
                    "PDCD1","FOXP3","HAVCR2","LAG3","CTLA4",
                    "CD7","FGFBP2","KLRF1","KLRB1","NCR1","FCGR3A", #NK
                    "CD14","CD68","CD163","S100A8","S100A9","MMP19", #MO
                    'CSF3R',
                    "APOC1","APOE","C1QA","C1QB", #Mac
                    "CLEC9A","CD1C", #DC
                    'TPSAB1','CPA3','CST3', # Mast
                    "ACTA2","COL1A2",'COL1A1',#Fibroblast
                    "MME",'FGF7',
                    "ENG","VWF","PECAM1" #Endo
                    )
DotPlot(sce,
        features = genes_to_check) + RotatedAxis()
sce[['cell_type']] = unname(cluster2celltype[sce@meta.data$seurat_clusters])

sce_epi <- subset(sce,cell_type == 'Epithelial_cells')
Epi_count <- GetAssayData(object = sce_epi,slot = 'count') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'SYMBOL')
Epi_count[1:5,1:5]
write.csv(Epi_count,file = 'for_CancerFinder.csv',row.names = F,col.names = T)
output <- fread('CancerFinder_output.csv',data.table = F)
sce_epi$cancerfinder_res <- output$predict
cells <- colnames(sce_epi)[sce_epi$group == 'tumor' &
                             sce_epi$cancerfinder_res == 1]
sce@meta.data[cells,'cell_type'] <- 'Malignant'

DimPlot(sce,reduction="umap",
        group.by="cell_type",
        pt.size=0.5)+
  ggtitle('Cell Type')




