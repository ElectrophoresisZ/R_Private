library(clustree)

load('LAM1_NORMAL_Reduction.rda')
LAM1_NORMAL_2 <- JoinLayers(LAM1_NORMAL)
LAM1_NORMAL_2 <- FindNeighbors(object = LAM1_NORMAL_2, dims = 1:30, reduction = 'pca')
LAM1_NORMAL_2 <- FindClusters(
  object = LAM1_NORMAL_2,
  resolution = c(seq(0.1,1.0,0.1))
)
clustree(LAM1_NORMAL_2@meta.data,prefix = 'RNA_snn_res.') +
  ggtitle("4a: Cluster Results by Different Resolution") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))

p1 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "RNA_snn_res.0.1",
              label = T,
              repel = F,
              shuffle = T) +  # 打乱顺序，避免遮盖
  ggtitle("Resolution = 0.1") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))
  
p2 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "RNA_snn_res.0.2",
              label = T,
              repel = F,
              shuffle = T) +
  ggtitle("Resolution = 0.2") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))

p3 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "RNA_snn_res.0.3",
              label = T,
              repel = F,
              shuffle = T) + 
  ggtitle("Resolution = 0.3") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))

p4 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "RNA_snn_res.0.4",
              label = T,
              repel = F,
              shuffle = T) +
  ggtitle("Resolution = 0.4") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))

p5 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "RNA_snn_res.0.5",
              label = T,
              repel = F,
              shuffle = T) +
  ggtitle("Resolution = 0.5") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))

p6 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "RNA_snn_res.1",
              label = T,
              repel = F,
              shuffle = T) +
  ggtitle("Resolution = 1") +
  theme(plot.title = element_text(face = 'bold',size = 16,hjust = 0.5))

p = p1+p2+p3+p4+p5+p6
p


LAM1_NORMAL_2$seurat_clusters <- LAM1_NORMAL_2$RNA_snn_res.0.1

UMAPPlot(LAM1_NORMAL_2,
         label = T,
         repel = F,
         group.by = 'seurat_clusters'
)+
  ggtitle('4c: UMAP Plot (Resolution = 0.1)')+
  theme_bw()+
  theme(
    plot.title = element_text(face = 'bold',hjust = 0.5),
    plot.background = element_rect(fill = 'transparent',color = NA)
  )+
  scale_color_igv()

Idents(LAM1_NORMAL_2) = 'seurat_clusters'
#save(LAM1_NORMAL_2, file = 'LAM1_NORMAL_Classified.rda')
