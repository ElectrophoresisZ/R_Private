
load('LAM1_NORMAL_Classified.rda')
load('ALL_MARKERS.rda')
library(readxl)
marker <- read_xlsx('Cell_marker_Human.xlsx')
markers <- marker$marker
top10 <- all.markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>% #以metadata中的‘cluster’来进行分类
  filter(p_val_adj < 0.05) %>%
  top_n(10,avg_log2FC) %>%
  group_by()


library(dplyr)
# 首先按 cluster 列升序排列
# 对于相同 cluster 的行，再按 avg_log2FC 列升序排列
plotdt <- all.markers %>%
  filter(gene %in% c(markers)) %>%
  arrange(cluster,avg_log2FC) %>%
  mutate(gene = factor(gene,levels = unique(gene),ordered = T))
plotdt_arrange <- plotdt %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup()

plotdt_up <- plotdt %>%
  filter(avg_log2FC > 0) %>%
  arrange(cluster,desc(avg_log2FC))
plotdt_down <- plotdt %>%
  filter(avg_log2FC < 0) %>%
  arrange(cluster,avg_log2FC)

top20 <- plotdt_up %>%
  group_by(cluster) %>% #以metadata中的‘cluster’来进行分类
  filter(p_val_adj < 0.05) %>%
  top_n(20,avg_log2FC) %>%
  group_by()

library(ggplot2)
ggplot(
  plotdt_arrange,
  aes(x = cluster,y = gene,
      size = avg_log2FC,
      color = pct.1))+
  geom_point()+
  geom_text(aes(label = cluster),size = 3,color = 'black')+
  theme_bw()+
  labs(title = "6a: Gene Expression Bubble Plot") +
  theme(
    plot.title = element_text(face = 'bold',hjust = 0.5,size = 16),
    plot.background = element_rect(fill = 'transparent',color = NA)
  )+
  scale_color_gradient2(low = 'olivedrab',high = 'salmon',
                        mid = 'yellow',midpoint = 0.5)
#----------------------------------------------------------------------------------

good_markers <- all.markers %>%
  filter(cluster == 0) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20) 
cellmarker <- good_markers$gene
cluster0 <- marker[which(marker$marker %in% cellmarker),]
table(cluster0$cell_name) #Myeloid cell 
# (Macrophage,Monocyte,Intermediate monocyte,Classical monocyte )

good_markers <- all.markers %>%
  filter(cluster == 1) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster1 <- marker[which(marker$marker %in% cellmarker),]
table(cluster1$cell_name) #Macrophage (probably M2)

good_markers <- all.markers %>%
  filter(cluster == 2) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster2 <- marker[which(marker$marker %in% cellmarker),]
table(cluster2$cell_name) #Epithelial cell 
# (Alveolar epithelial cell Type 2,Alveolar macrophage)

good_markers <- all.markers %>%
  filter(cluster == 3) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster3 <- marker[which(marker$marker %in% cellmarker),]
table(cluster3$cell_name) #T cell ;Natural killer cell
# (Natural killer T (NKT) cell, Regulatory T(Treg) cell,Cytotoxic T cell,CD8+ T cell,CD4+ T cell)

good_markers <- all.markers %>%
  filter(cluster == 4) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster4 <- marker[which(marker$marker %in% cellmarker),]
table(cluster4$cell_name) #Endothelial cell

good_markers <- all.markers %>%
  filter(cluster == 5) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster5 <- marker[which(marker$marker %in% cellmarker),]
table(cluster5$cell_name) #Endothelial cell

good_markers <- all.markers %>%
  filter(cluster == 6) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster6 <- marker[which(marker$marker %in% cellmarker),]
table(cluster6$cell_name) #Fibroblast 
#(probably include Stromal cell, Smooth muscle cell,Pericyte,Myofibroblast,Mesenchymal cell)

good_markers <- all.markers %>%
  filter(cluster == 7) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster7 <- marker[which(marker$marker %in% cellmarker),]
table(cluster7$cell_name) #Endothelial cell
#(probably Lymphatic endothelial cell, include Fibroblast)

good_markers <- all.markers %>%
  filter(cluster == 8) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster8 <- marker[which(marker$marker %in% cellmarker),]
table(cluster8$cell_name) #Ciliated epithelial cell

good_markers <- all.markers %>%
  filter(cluster == 9) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster9 <- marker[which(marker$marker %in% cellmarker),]
table(cluster9$cell_name) #?

good_markers <- all.markers %>%
  filter(cluster == 10) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster10 <- marker[which(marker$marker %in% cellmarker),]
table(cluster10$cell_name) #Red blood cell (erythrocyte)

good_markers <- all.markers %>%
  filter(cluster == 11) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster11 <- marker[which(marker$marker %in% cellmarker),]
table(cluster11$cell_name) #Neural progenitor cell;MKI67+ progenitor cell

good_markers <- all.markers %>%
  filter(cluster == 12) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster12 <- marker[which(marker$marker %in% cellmarker),]
table(cluster12$cell_name) #Macrophage 

good_markers <- all.markers %>%
  filter(cluster == 13) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster13 <- marker[which(marker$marker %in% cellmarker),]
table(cluster13$cell_name) #T cell

good_markers <- all.markers %>%
  filter(cluster == 14) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster14 <- marker[which(marker$marker %in% cellmarker),]
table(cluster14$cell_name) #Plasma cell
# (Memory B cell,Plasmablast)

LAM1_NORMAL_2

LAM1_NORMAL_2$cell_type = case_when(
  LAM1_NORMAL_2$seurat_clusters %in% c(2,8) ~ 'Epithelial cell',
  LAM1_NORMAL_2$seurat_clusters %in% c(9) ~ 'Unknown',
  LAM1_NORMAL_2$seurat_clusters %in% c(4,5,7) ~ 'Endothelial cell',
  LAM1_NORMAL_2$seurat_clusters %in% c(11) ~ 'Neural progenitor cell',
  LAM1_NORMAL_2$seurat_clusters %in% c(6) ~ 'Fibroblast ',
  LAM1_NORMAL_2$seurat_clusters %in% c(10) ~ 'Red blood cell',
  LAM1_NORMAL_2$seurat_clusters %in% c(9) ~ 'Unknown',
  LAM1_NORMAL_2$seurat_clusters %in% c(1,12) ~ 'Macrophage',
  LAM1_NORMAL_2$seurat_clusters %in% c(3,13) ~ 'T cells',
  LAM1_NORMAL_2$seurat_clusters %in% c(14) ~ 'Plasma cell',
  LAM1_NORMAL_2$seurat_clusters %in% c(0) ~ 'Myeloid cell',
  TRUE ~ 'Unknown'
)

#细胞注释
p1 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "seurat_clusters",
              label = T,
              repel = F,
              shuffle = T) + labs(title = "6b: Annotation by Original Clusters") +
  theme(plot.title = element_text(face='bold',hjust = 0.5,size = 16))
p2 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "cell_type",
              label = T,
              repel = F,
              shuffle = T) + labs(title = "6c: Annotation by Integrated Cell Type") +
  theme(plot.title = element_text(face='bold',hjust = 0.5,size = 16))
p3 <- DimPlot(LAM1_NORMAL_2,reduction='umap',group.by = "type",
              label = T,
              repel = F,
              shuffle = T) +labs(title = "6d: Annotation by Cell Type") +
  theme(plot.title = element_text(face='bold',hjust = 0.5,size = 16))
p = p1+p2+p3
p

plotdata = LAM1_NORMAL_2@meta.data

#分组细胞数量分布图
ggplot(
  plotdata,
  aes(x=cell_type,fill = type)
)+
  geom_bar(position = 'dodge')+  #并排显示分组条形图
  scale_fill_igv(alpha = 0.7) +
  labs(title = "6e: Cell Count of Different Types",x="") +
  theme(plot.title = element_text(face='bold',hjust = 0.5,size = 21),
        axis.title = element_text(face='bold',hjust = 0.5,size = 16),
        axis.text.y= element_text(face='bold',hjust = 0.5,size = 12)) +
  coord_flip()

# 百分比堆叠图
ggplot(
  plotdata,
  aes(x=type,fill = cell_type)
)+
  geom_bar(position = 'fill')+ #百分比堆叠
  scale_fill_igv(alpha = 0.7)+
  labs(title = "6f: Cell Percentage of Different Types",y="Percentage") +
  theme(plot.title = element_text(face='bold',hjust = 0.5,size = 21),
        axis.title = element_text(face='bold',hjust = 0.5,size = 16),
        axis.text.x= element_text(face='bold',hjust = 0.5,size = 12))

UMAPPlot(LAM1_NORMAL_2,
         label = T,
         repel = T,
         split.by = 'type',  #按指定条件分割图形
         group.by = 'cell_type'
)+
  ggtitle('6g: UMAP by Type')+
  theme_bw()+
  theme(
    plot.title = element_text(face = 'bold',hjust = 0.5),
    plot.background = element_rect(fill = 'transparent',color = NA)
  )+
  scale_color_igv()
#save(LAM1_NORMAL_2, file = "LAM1_NORMAL_Annotated.rda")
