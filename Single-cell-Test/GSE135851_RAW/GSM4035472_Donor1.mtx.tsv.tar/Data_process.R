library(Seurat)
library(patchwork)
normal.data <- Read10X('Normal_Lung/')
dim(normal.data)
normal <- CreateSeuratObject(counts = normal.data,
                           project = "normal_lung",
                           min.cells = 3,   #只保留那些在至少 3 个细胞中表达的基因
                           min.features = 200) #只保留那些表达至少 200 个基因的细胞
normal
normal[['percent.mt']] <- PercentageFeatureSet(normal,pattern = '^MT-')
VlnPlot(normal,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

# ncount:细胞中所有基因的表达量（原始计数）的总和
# nfeature:表达水平大于 0 的基因数
plot1 <- FeatureScatter(normal,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(normal,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2

normal <- subset(normal, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
normal <- NormalizeData(normal)
normal
normal <- FindVariableFeatures(normal,selection.method = 'vst',
                               nfeatures = 2000)
top10 <- head(VariableFeatures(normal),10)
plot3 <- VariableFeaturePlot(normal)
plot4 <- LabelPoints(plot3, points = top10, repel = T)
plot4

all.genes <- rownames(normal)
normal <- ScaleData(normal, features = all.genes)

#对高变基因降维
normal <- RunPCA(normal,features = VariableFeatures(normal))
print(normal[['pca']],dims = 1:5,nfeatures = 5) #前五个主成分和前五个基因
VizDimLoadings(normal,dims = 1:2,reduction = 'pca') #指定维度的基因加载情况
DimPlot(normal,reduction = 'pca')
DimHeatmap(normal,dims = 1:2, 
           cells = 500, # 随机抽样 500 个细胞
           balanced = T) # 同时显示正负 loadings 最高的基因

# JackStraw 置换检验
# 评估主成分分析中每个主成分的统计显著性
normal <-JackStraw(normal,num.replicate = 100)
# 计算每个主成分的总体显著性 p 值
normal <- ScoreJackStraw(normal,dims = 1:20)
JackStrawPlot(normal,dims = 1:15)

#KNN聚类
normal <- FindNeighbors(normal, dims = 1:10)
normal <- FindClusters(normal, resolution = 0.5) #控制着分群的“粒度”
head(Idents(normal))
#非线性降维方法 (UMAP/tSNE)分群
normal <- RunUMAP(normal,dims = 1:10)
DimPlot(normal, reduction = "umap",label = T)

all.markers <- FindAllMarkers(object = normal, 
                              only.pos = TRUE, #只返回上调的差异表达基因
                              min.pct = 0.25, #至少 25% 的细胞表达的基因才会被考虑为差异表达基因
                              test.use = "wilcox",
                              logfc.threshold = 0.25)
#save(all.markers,file = 'normal_markers.Rdata')

library(dplyr)
top2gene <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) # 选取avglogFC最大的基因
# 基因分布图
FeaturePlot(normal, 
            features = top2gene$gene[which(top2gene$avg_log2FC>5)])

top10gene <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10,  order_by= avg_log2FC)
#基因在不同簇中的热图
DoHeatmap(normal, features = top10gene$gene) + NoLegend()
top10gene_select <- top10gene[which(top10gene$avg_log2FC>5),]
DotPlot(normal, features = unique(top10gene_select$gene))+RotatedAxis()

saveRDS(normal,file = 'normal_lung_final.rds')


