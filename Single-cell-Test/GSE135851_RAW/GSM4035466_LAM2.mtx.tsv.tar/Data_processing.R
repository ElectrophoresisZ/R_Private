library(Seurat)
library(patchwork)
LAM2.data <- Read10X('LAM2/')
dim(LAM2.data)
LAM2 <- CreateSeuratObject(counts = LAM2.data,
                           project = 'LAM',
                           min.cells = 3,
                           min.features = 200)
LAM2
LAM2[['percent.mt']] <- PercentageFeatureSet(LAM2,pattern = '^MT-')
#QC_violin_plot
VlnPlot(LAM2, features = colnames(LAM2@meta.data)[2:4], ncol = 3)
colnames(LAM2@meta.data)[2:4]

#Feature_scatter_plot
plot1 <- FeatureScatter(LAM2,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(LAM2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1+plot2

LAM2 <- subset(LAM2, subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<5)
LAM2 <- NormalizeData(LAM2)
LAM2
LAM2 <- FindVariableFeatures(LAM2,selection.method = 'vst',nfeatures = 2000)
top10 <- head(VariableFeatures(LAM2),10)
top10
#Variable_feature_plot
plot3 <- VariableFeaturePlot(LAM2)
plot4 <- LabelPoints(plot3,points = top10,repel = T)
plot4

all.genes <- rownames(LAM2)
LAM2 <- ScaleData(LAM2,features = all.genes)

LAM2 <- RunPCA(LAM2,features = VariableFeatures(LAM2))
#PCA_loading_plot
VizDimLoadings(LAM2,dims = 1:4,reduction = 'pca')
#PCA_scatter_plot
DimPlot(LAM2,reduction = 'pca')
#PCA_heatmap
DimHeatmap(LAM2,dims = 1:2,
           cells = 500,
           balanced = T)
LAM2 <- JackStraw(LAM2,reduction = 'pca',num.replicate = 100)
LAM2 <- ScoreJackStraw(LAM2,dims = 1:20)
#Jackstraw_plot
JackStrawPlot(LAM2,dims = 1:15)

LAM2 <- FindNeighbors(LAM2,dims = 1:10)
LAM2 <- FindClusters(LAM2,resolution = 0.5)
head(Idents(LAM2))
LAM2 <- RunUMAP(LAM2,dims=1:10)
#UMAP_plot
DimPlot(LAM2,reduction = 'umap',label = T)

all.markers <- FindAllMarkers(object = LAM2,
                              only.pos = T,
                              min.pct = 0.25,
                              test.use = 'wilcox',
                              logfc.threshold = 0.25)
#save(all.markers,file = 'all.markers.Rdata')
library(dplyr)
top2gene <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n=2,order_by = avg_log2FC)
#Gene_exp_UMAP_plot
FeaturePlot(LAM2,features = top2gene$gene[which(as.numeric(top2gene$avg_log2FC)>6)])

top10gene <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n=10,order_by = avg_log2FC)
#top10gene <- top10gene[order(as.numeric(top10gene$avg_log2FC),decreasing = T),]
#Gene_dotplot
DotPlot(LAM2,features = unique(top10gene$gene)[1:20])+RotatedAxis()

library(SingleR)
library(celldex)
library(scrapper)
ref <- HumanPrimaryCellAtlasData()
Idents(LAM2)
LAM2_singler <- GetAssayData(LAM2,slot = 'data')
clusters <- LAM2$seurat_clusters
pred.annotation <- SingleR(test = LAM2_singler,
                           ref = ref,
                           labels = ref$label.fine,
                           method = 'cluster',
                           clusters = clusters)
table(pred.annotation$labels)
#SingleR_score_heatmap
plotScoreHeatmap(pred.annotation,clusters = pred.annotation$labels)
newcluster <- pred.annotation$labels
names(newcluster) <-levels(LAM2)
LAM2 <- RenameIdents(LAM2, newcluster)
#Annotated_UMAP_plot
DimPlot(LAM2,reduction = 'umap',label = T,repel = T)
