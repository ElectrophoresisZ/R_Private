load("LAM1_NORMAL_Annotated.rda")
Idents(LAM1_NORMAL_2) <- "cell_type"
cell.markers <- FindAllMarkers(
  object = LAM1_NORMAL_2,
  only.pos = F, #保留上调与下调的基因
  test.use = 'wilcox', #检验方法是wilcox检验
  slot = 'data', #默认使用data而不是counts
  min.pct = 0.25, #在至少25%的细胞内表达
  logfc.threshold = 0.25 #设置log2FC的阈值
)
#save(cell.markers,file = "MANUAL_ALL_MARKERS.rda")


degs = FindMarkers(
  LAM1_NORMAL_2,
  logfc.threshold = 0.5,
  only.pos = F,
  ident.1 = '',
  ident.2 = '',#默认使用ident1 vs ident2，而不是相反
  group.by = 'cell_type'
) %>%
  mutate(gene=rownames(.))

Epithelial_cell_subtype = subset(LAM1_NORMAL_2 , cell_type == 'Epithelial cell')
Epithelial_cell_degs = FindMarkers(
  Epithelial_cell_subtype,
  logfc.threshold = 0.5,
  only.pos = F,
  ident.1 = 'LAM',
  ident.2 = 'normallung',
  group.by = 'type'
) %>%
  mutate(gene=rownames(.))

change = as.factor(ifelse(Epithelial_cell_degs$p_val_adj < 0.05 & abs(Epithelial_cell_degs$avg_log2FC) > 1,
                          ifelse(Epithelial_cell_degs$avg_log2FC > 1 ,'Up','Down'),'No change'))

filtered_Epithelial_cell_degs = Epithelial_cell_degs %>%
  filter(p_val_adj<1e-100)%>%
  filter(abs(avg_log2FC)>1)
Epithelial_cell_degs$label <- ifelse(Epithelial_cell_degs$gene %in% filtered_Epithelial_cell_degs$gene, Epithelial_cell_degs$gene, NA)

ggplot(Epithelial_cell_degs, aes(avg_log2FC, -log10(p_val_adj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+  #设置基础字体大小
  labs(title = "7a: Volcano Plot of Epithelial Cells") +
  geom_text_repel(aes(label = label,fontface = 'bold'), 
                  size = 3,
                  box.padding = 0.2,  # 标签与点之间的间距
                  max.overlaps = 10) + # 超过10个重叠时只显示部分标签
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No change" = "gray"))+
  theme(panel.grid = element_blank(), #移除面板网格线
        legend.position = 'right',
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 16))+
  scale_x_continuous(limits = c(-max(abs(Epithelial_cell_degs$avg_log2FC)), max(abs(Epithelial_cell_degs$avg_log2FC)))) #设置x轴范围对称

celldiff_top3 <- cell.markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 3) %>%
  slice_max(avg_log2FC, n = 3) %>%  # 选择logFC最高的3个基因
  group_by()

VlnPlot(LAM1_NORMAL_2,features = celldiff_top3$gene,
        group.by = "cell_type",
        split.plot = T, #创建分割图
        split.by = 'type',
        cols = c('olivedrab','salmon'),
        stack = T, #创建分割图
        flip = T)  +#将x轴和y轴翻转
  labs(title = "7b: Violin Plot of Annotated Cells") +
  theme(plot.title = element_text(face = 'bold',hjust = 0.5,size = 16))

