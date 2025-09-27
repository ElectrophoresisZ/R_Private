load('LAM_NORMAL_filtered.rda')
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)
LAM_NORMAL <- LAM1_NORMAL_scaled
LAM_NORMAL <- NormalizeData(object = LAM_NORMAL)
LAM_NORMAL <- FindVariableFeatures(object = LAM_NORMAL)
VariableFeaturePlot(LAM_NORMAL,cols=c('grey','steelblue'))+ #高变基因的点颜色:蓝色
  scale_y_log10() + #将y轴转换为对数尺度
  labs(title = '3a: Variable Feature Plot') 

par(mfrow=c(1,1))
LAM_NORMAL <- ScaleData(LAM_NORMAL)
hist(colSums(LAM_NORMAL,slot='scale.data'),
     breaks = 50,
     main = '3b: Scaled Gene Expression Distribution',
     xlab = "Total scaled expression per cell",    # x轴：每个细胞的总标准化表达量
     ylab = "Number of cells")

LAM_NORMAL <- RunPCA(LAM_NORMAL)
ElbowPlot(LAM_NORMAL,ndims = 30) +
  labs(title = "3c: Principal component elbow diagram")

LAM_NORMAL <- RunTSNE(LAM_NORMAL,dims = 1:30,reduction = 'pca',reduction.name = 'tsne')
LAM_NORMAL <- RunUMAP(object = LAM_NORMAL, dims = 1:30, reduction = 'pca',reduction.name = 'umap')
LAM_NORMAL <- FindNeighbors(object = LAM_NORMAL, dims = 1:30, reduction = 'pca')#根据前30个主成分构建细胞间邻接图
LAM_NORMAL <- FindClusters(object = LAM_NORMAL,resolution = 0.2, cluster.name = 'unintegrated_clusters')#进行聚类（预聚类）

DimPlot(object = LAM_NORMAL, reduction = "pca") +
  labs(title = "3d: PCA for Dimensionality Reduction")
DimPlot(object = LAM_NORMAL, reduction = "tsne") +
  labs(title="3e: Tsne for Dimensionality Reduction")
DimPlot(object = LAM_NORMAL, reduction = "umap") +
  labs(title="3f: Umap for Dimensionality Reduction")

library(patchwork)
p1 <- DimPlot(object = LAM_NORMAL, reduction = 'umap',
        group.by = 'type',
        shuffle = T) + labs(title = "3g: Classification by Cell Type")

p2 <- DimPlot(object = LAM_NORMAL, reduction = 'umap',
        group.by = 'unintegrated_clusters',
        shuffle = T) + labs(title = "3h: Classification by Unintegrated Clusters")
plot_grid(p1,p2,nrow = 1,ncol = 2)
LAM1_NORMAL <- LAM_NORMAL
#save(LAM1_NORMAL,file = 'LAM1_NORMAL_Reduction.rda')
