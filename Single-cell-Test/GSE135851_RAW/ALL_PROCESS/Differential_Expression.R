
load('LAM1_NORMAL_Classified.rda')
all.markers <- FindAllMarkers(LAM1_NORMAL_2,
                              only.pos = F,  #保留上调与下调的基因
                              test.use = 'wilcox',
                              slot = 'data',     #默认使用data而不是counts
                              min.pct = 0.25,   #该基因在至少25%的细胞内表达
                              logfc.threshold = 0.25 
                              )
#save(all.markers,file = 'ALL_MARKERS.rda')

top3 <- all.markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>% #以metadata中的‘cluster’来进行分类
  filter(p_val_adj < 0.05) %>%
  top_n(3,avg_log2FC) %>%
  group_by()

DoHeatmap(
  object = LAM1_NORMAL_2,
  features = top3$gene,
  label = TRUE
)+
  scale_fill_gradient(low = 'lightgrey',high = 'salmon') +
  labs(title = "5a: Gene Heatmap by Different Clusters") +
  theme(plot.title = element_text(face="bold",size=16,hjust = 0.5))
  
DotPlot(
  object = LAM1_NORMAL_2,
  group.by = 'seurat_clusters',
  cols = c('lightgrey','salmon'),
  features = unique(top3$gene)
)+
  coord_flip() +
  labs(title = "5b: Gene Dotplot by Different Clusters") +
  theme(plot.title = element_text(face="bold",size=16,hjust = 0.5))

VlnPlot(
  object = LAM1_NORMAL_2,
  group.by = 'seurat_clusters',
  #fill.by = 'feature',
  stack = T,
  flip = T,
  features = unique(top3$gene)
)+
  scale_fill_igv() +
  labs(title = "5c: Gene Violinplot by Different Clusters") +
  theme(plot.title = element_text(face="bold",size=16,hjust = 0.5))


FeaturePlot(
  LAM1_NORMAL_2,
  features = top3$gene[37:45],#可以选择感兴趣的基因进行绘图
  order = T,
  alpha = 0.3,
  reduction = 'umap'
) 

library(patchwork)
p1 <- DimPlot(object = LAM1_NORMAL_2, reduction = 'umap',
              group.by = 'type',
              shuffle = T) + labs(title = "5e: Classification by Cell Type")

p2 <- DimPlot(object = LAM1_NORMAL_2, reduction = 'umap',
              group.by = 'seurat_clusters',
              shuffle = T) + labs(title = "5f: Classification by Seurat Clusters")
plot_grid(p1,p2,nrow = 1,ncol = 2)


