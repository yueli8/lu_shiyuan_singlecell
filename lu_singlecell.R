#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
library(dplyr)
library(patchwork)
library(Seurat)
library(plyr)
dir="~/lu_shiyuan/filtered_feature_bc_matrix"
list.files(dir)
pbmc.data<-Read10X(data.dir=dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "~/lu_shiyuan/lu_cluster.rds")

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)
cluster6.markers <- FindMarkers(pbmc, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 5)
cluster7.markers <- FindMarkers(pbmc, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 5)
cluster8.markers <- FindMarkers(pbmc, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 5)
cluster9.markers <- FindMarkers(pbmc, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 5)

new.cluster.ids <- c("Natural Killer T Cell", "Natural Killer T Cell", "Natural Killer T Cell", "Natural Killer T Cell", "Progenitor Cell",
                     "Progenitor Cell","Germ Cell","Natural Killer T Cell","Natural Killer T Cell","T Helper2 Cell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "~/lu_shiyuan/lu_cluster_id.rds")

setwd("/home/hp/lu_shiyuan")
lu_cluster_id<-readRDS(file="lu_cluster_id.rds")
pbmc<-readRDS(file="lu_cluster_id.rds")

Natural_Killer <-subset(lu_cluster_id, idents=c('Natural Killer T Cell'))
DimPlot(Natural_Killer, reduction = "umap")
saveRDS(Natural_Killer, file="Natural_Killer.rds")

T_Helper2 <-subset(lu_cluster_id, idents=c('T Helper2 Cell'))
DimPlot(T_Helper2, reduction = "umap")
saveRDS(T_Helper2, file="T_Helper2.rds")

Germ_Cell <-subset(lu_cluster_id, idents=c('Germ Cell'))
DimPlot(Germ_Cell, reduction = "umap")
saveRDS(Germ_Cell, file="Germ_Cell.rds")

Progenitor_Cell <-subset(lu_cluster_id, idents=c('Progenitor Cell'))
DimPlot(Progenitor_Cell, reduction = "umap")
saveRDS(Progenitor_Cell, file="Progenitor_Cell.rds")

Natural_Killer<-readRDS(file="Natural_Killer.rds")

#a<-DoHeatmap(subset(Natural_Killer,downsample=50000),size=5)
#a1<-a$data
#write.table(a1,"a1")
#a2<-a1[c('Feature','Expression')]
#na.rm=TURE为忽略缺失数据的意思，sum(！is.na())为计算非缺失数据的频数。
#a2_average<-ddply(a2,c("Feature"),expression=mean(Expression,na.rm=TRUE))
#a1_unique<-unique(a2_average)
#write.table(a1_unique,"a1_unique")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

VlnPlot(pbmc, features = c("CD52","CD96","CDC25B","EEF2","HLA-A",
                           "HLA-B", "HLA-DRA", "HLA-DRB1", "HLA-DRB5","IL32",
                           "LTB","PTPRC","RPL39","S100A6","SH3BGRL3","UCP2"))


NKT.markers <- FindMarkers(pbmc, ident.1 = "Natural Killer T Cell", min.pct = 0.25)
head(NKT.markers, n = 12)

NKT.markers <- FindMarkers(pbmc, ident.1 = "Natural Killer T Cell", min.pct = 0.25)
head(NKT.markers, n = 1000)

write.table(NKT.markers,file="NKT.markers")

FeaturePlot(pbmc, features = c("CD52","CD96","CDC25B", "EEF2",  
                               "HLA-A", "HLA-B", "HLA-DRA", "HLA-DRB1", "HLA-DRB5","IL32",
                               "LTB","PTPRC","RPL39","S100A6","SH3BGRL3","UCP2"))




