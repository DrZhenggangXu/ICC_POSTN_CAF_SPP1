getwd()
rm(list = ls())
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
library(CellChat)
data.input = sce@assays$RNA@data
meta.data =  sce@meta.data
unique(meta.data$cell_type)
data.input <- data.input[,rownames(meta.data)]
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = meta.data)
levels(cellchat@idents)
cellchat <- setIdent(cellchat, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,
                              raw.use = FALSE,
                              population.size = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights")
netVisual_bubble(cellchat,
                 vjust.x = 0,
                 font.size = 16,
                 font.size.title = 16,
                 show.legend = F
                 )

