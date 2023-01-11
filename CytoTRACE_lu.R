#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("slingshot")

#library("devtools")
#install_github("jw156605/SLICER")

#install.packages("devtools")
#devtools::install_github("dynverse/dyno")

setwd("~/lu_shiyuan")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="lu_cluster_id.rds")
DimPlot(seurat, label = T) + NoLegend()
table(Idents(seurat))

##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)

result01<-CytoTRACE(a,ncores=4)#ncores only can be 4 , can not be 8.

plotCytoGenes(result01, numOfGenes = 10)#first figure

cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')
p1

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('int.umap')

p1|p2

p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))

colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="Natural Killer T",
                                                 "1"="Natural Killer T",
                                                 "2"="Natural Killer T",
                                                 "3"="Natural Killer T",
                                                 "4"="Progenitor",
                                                 "5"="Progenitor",
                                                 "6"="Germ Cell",
                                                 "7"="Natural Killer T",
                                                 "8"="Natural Killer T",
                                                 "9"="T Helper2 ")


colnames(colData(cds))
head(colData(cds))

a<-colData(cds)
write.csv(a,"a")

table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

plotCytoTRACE(result01, phenotype = phe)

plotCytoTRACE(result01, phenotype = phe,gene = "CD52",outputDir = "CD52")
plotCytoTRACE(result01, phenotype = phe,gene = "CD96",outputDir = "CD96")
plotCytoTRACE(result01, phenotype = phe,gene = "CDC25B",outputDir = "CDC25B")
plotCytoTRACE(result01, phenotype = phe,gene = "EEF2",outputDir = "EEF2")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-A",outputDir = "HLA-A")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-B",outputDir = "HLA-B")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-DRA",outputDir = "HLA-DRA")
plotCytoTRACE(result01, phenotype = phe,gene = "HLA-DRB1",outputDir = "HLA-DRB1")
plotCytoTRACE(result01, phenotype = phe,gene = "IL32",outputDir = "IL32")
plotCytoTRACE(result01, phenotype = phe,gene = "LTB",outputDir = "LTB")
plotCytoTRACE(result01, phenotype = phe,gene = "PTPRC",outputDir = "PTPRC")
plotCytoTRACE(result01, phenotype = phe,gene = "RPL39",outputDir = "RPL39")
plotCytoTRACE(result01, phenotype = phe,gene = "S100A6",outputDir = "S100A6")
plotCytoTRACE(result01, phenotype = phe,gene = "SH3BGRL3",outputDir = "SH3BGRL3")
plotCytoTRACE(result01, phenotype = phe,gene = "UCP2",outputDir = "UCP2")
