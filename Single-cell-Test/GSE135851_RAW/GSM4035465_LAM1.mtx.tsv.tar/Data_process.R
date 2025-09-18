library(Seurat)
library(patchwork)
pbmc.data <- Read10X('LAM1/')
dim(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.data,
                           project = "pbmck",
                           min.cells = 3,
                           min.features = 200)
pbmc
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc,pattern = '^MT-')

# orig.ident: 项目标识符（设置的 project = "pbmc3k"）
# nCount_RNA: 每个细胞的总UMI数（所有基因的计数之和）
# nFeature_RNA: 每个细胞中检测到的非零基因数量
VlnPlot(pbmc,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3)
plot1 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
pbmc

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = 'vst', # 方差稳定变换
                             nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc),10)
top10

plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3,points = top10,repel = T) # 给前10个加名字
plot3 + plot4


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) 


#PCA降维
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))
print(pbmc[['pca']], dims = 1:5,nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#KNN聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) #控制着分群的“粒度”
head(Idents(pbmc))


#非线性降维方法 (UMAP/tSNE)分群
pbmc <- RunUMAP(pbmc,dims = 1:10)
DimPlot(pbmc, reduction = "umap",label = T)

#saveRDS(pbmc,file = 'pbmc_LAM1.rds')

all.markers <- FindAllMarkers(object = pbmc, 
                              only.pos = TRUE, 
                              min.pct = 0.25,
                              test.use = "wilcox",
                              logfc.threshold = 0.25)
# 查看每个类群中前2个marker基因
library(dplyr)
top2gene <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
VlnPlot(pbmc,features = c("RP11-598F7.3","SLC6A14","IGLL5"))
top2gene$gene
FeaturePlot(pbmc, 
            features = top2gene$gene[which(top2gene$pct.1>0.5)])

top10 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top20 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

# 展示每个类群前10个基因
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
top10_updata <- top10[order(top10$avg_log2FC,decreasing = T),]
DotPlot(pbmc, features = unique(top10$gene)[1:20])+RotatedAxis()

library(readxl)
cellmarker <- read_xlsx('Cell_marker_Human.xlsx')
overlap <- inner_join(top10, cellmarker, by = c("gene" = "marker"))
overlap_sim <- data.frame(gene=overlap$gene,cluster=overlap$cluster,cell=overlap$cell_name)
new.cluster.ids <- c("Alveolar Epithelial Cell Type II ", 'Macrophage','','Neutrophil',
                     'Club cell','T cell','M2 macrophage','Lymphatic Endothelial Cell',
                     'Monocyte','Plasma cell','B cell','M2 Macrophage','Endothelial cell') 

names(new.cluster.ids) <- levels(pbmc)

pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
#saveRDS(pbmc,file = 'pbmc_LAM1_final.rds')
pbmc <- readRDS('pbmc_LAM1_final.rds')
