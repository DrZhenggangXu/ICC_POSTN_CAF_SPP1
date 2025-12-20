getwd()
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SingleR)
library(devtools)
library(harmony)
library(DoubletFinder)
library(multtest)
library(Rcpp)
library(ggplot2)
library(clustree)
library(org.Hs.eg.db)
library(clusterProfiler)
library(RColorBrewer)
require(monocle)
sce <- readRDS('sce_fibro.rds')
Counts <- GetAssayData(object = sce,slot = 'counts',assay="RNA")
gene_annotation <- data.frame(gene_short_name = row.names(Counts),
                              biotype=rep("protein_coding",nrow(Counts)))
rownames(gene_annotation)<-rownames(Counts)
sample_sheet<-as.data.frame(sce@meta.data)
sample_sheet<-cbind(rownames(sample_sheet),sample_sheet)
colnames(sample_sheet)[1]<-"cell"
sample_sheet<-as.data.frame(t(subset(t(sample_sheet),select=colnames(Counts))))
pd <- new('AnnotatedDataFrame', data = sample_sheet) 
fd <- new('AnnotatedDataFrame', data = gene_annotation)
cds <- newCellDataSet(cellData = Counts,  
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 1.0)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed > nrow(sample_sheet) * 0.01))
diff <-differentialGeneTest(cds[expressed_genes,],
                            fullModelFormulaStr="~cell_type",
                            reducedModelFormulaStr = "~orig.ident",
                            relative_expr=TRUE,cores=4)
ordering_genes <- row.names (subset(diff, qval < 0.001))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds,
                       max_components = 2,
                       reduction_method = 'DDRTree',
                       residualModelFormulaStr = "~orig.ident")
cds <- orderCells(cds, reverse = T)
plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE,cell_size = 1.0) +
  theme(legend.position="right")
plot_cell_trajectory(cds,color_by="cell_type",size=1,show_backbone=TRUE,cell_size = 1.0) +
  theme(legend.position="right")
plot_cell_trajectory(cds, color_by = "State", size=1,show_backbone=TRUE)
plot_genes_in_pseudotime(cds[row.names(fData(cds)) == "POSTN", ],
                         color_by = "cell_type") +
  theme(legend.position="right")

