library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)

##数据读取----------------------------------------------------------------------
LAM1.data <- Read10X('LAM1/')
LAM2.data <- Read10X('LAM2/')
LAM3.data <- Read10X('LAM3/')
par(mfrow=c(1,3))
hist(log2(LAM1.data@x),breaks=100)
hist(log2(LAM2.data@x),breaks=100)
hist(log2(LAM3.data@x),breaks=100)
par(mfrow=c(1,1))

LAM1 <- CreateSeuratObject(counts = LAM1.data,
                           project = "LAM1",
                           min.cells = 3,
                           min.features = 200)
LAM2 <- CreateSeuratObject(counts = LAM2.data,
                           project = "LAM2",
                           min.cells = 3,
                           min.features = 200)
LAM3 <- CreateSeuratObject(counts = LAM3.data,
                           project = "LAM3",
                           min.cells = 3,
                           min.features = 200)
LAM1;LAM2;LAM3
LAM12 <- merge(LAM1,LAM2)
LAM123 <- merge(LAM12,LAM3)

##质控---------------------------------------------------------------------------
LAM123$Percent.Mito <- PercentageFeatureSet(LAM123,pattern = '^MT-')
LAM123$Percent.ERCC <- PercentageFeatureSet(LAM123,pattern = '^ERCC')
LAM123$Percent.Ribo <- PercentageFeatureSet(LAM123,pattern = '^RP[SL]')

par(mfrow=c(1,3))
hist(LAM123$Percent.Mito)
hist(LAM123$Percent.ERCC)
hist(LAM123$Percent.Ribo)

p1 <- VlnPlot(LAM123,features = 'nFeature_RNA',group.by = 'orig.ident',
              alpha = 0.1,pt.size = 0.01) + # 控制每个数据点的大小
  geom_hline(yintercept = 2000,color = 'red')+  # 添加一条水平线;指定水平线的y轴位置为2000
  scale_fill_igv() # 提供了一套预设的颜色方案

p2 <- VlnPlot(LAM123,features = 'nCount_RNA',group.by = 'orig.ident',
              alpha = 0.1,pt.size = 0.01)+
  geom_hline(yintercept = 5000,color = 'red')+
  scale_fill_igv()

p3 <- VlnPlot(LAM123,features = 'Percent.Ribo',group.by = 'orig.ident',
              alpha = 0.1,pt.size = 0.01)+
  geom_hline(yintercept = 30,color = 'red')+
  scale_fill_igv()

p4 <- VlnPlot(LAM123,features = 'Percent.Mito',group.by = 'orig.ident',
              alpha = 0.1,pt.size = 0.01)+
  geom_hline(yintercept = 10,color = 'red')+
  scale_fill_igv()

library(cowplot)
plot_grid(p1, p2, p3, p4,
          nrow = 2, ncol = 2)

correlation <- cor(LAM123$nFeature_RNA, LAM123$nCount_RNA)
cor_text <- paste0("R = ", round(correlation, 3))

FeatureScatter(object = LAM123,
               group.by = 'orig.ident',
               raster = F, #是否将点栅格化
               shuffle = T, #是否随机打乱点的绘制顺序
               pt.size = 0.05,
               feature1 = "nFeature_RNA",
               feature2 = "nCount_RNA")+
  scale_color_igv()+
  guides(color = guide_legend(override.aes = list(size = 4))) + #让图例中的点变大
  labs(title = paste0("Correlation Coefficient: ", cor_text)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

correlation <- cor(LAM123$Percent.Ribo, LAM123$Percent.Mito)
cor_text <- paste0("R = ", round(correlation, 3))

FeatureScatter(object = LAM123,
               group.by = 'orig.ident',
               raster = F, #是否将点栅格化
               shuffle = T, #是否随机打乱点的绘制顺序
               pt.size = 0.05,
               feature1 = "Percent.Ribo",
               feature2 = "Percent.Mito")+
  scale_color_igv()+
  guides(color = guide_legend(override.aes = list(size = 4))) + #让图例中的点变大
  labs(title = paste0("Correlation Coefficient: ", cor_text)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))


LAM123_filter <- subset(LAM123,nCount_RNA>1000)
LAM123_filter <- subset(LAM123_filter,nFeature_RNA>500)
LAM123_filter <- subset(LAM123_filter,Percent.Mito<10)
LAM123_filter <- subset(LAM123_filter,Percent.Ribo<40)
LAM123_filter
#save(LAM123_filter,file = 'LAM123_filter.rda')

##降维聚类-----------------------------------------------------------------------
LAM123_Integrate <- LAM123_filter |> 
  Seurat::NormalizeData() |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
  ScaleData() |>
  RunPCA()

library(harmony)
LAM123_Integrate <- RunHarmony(LAM123_Integrate,group.by.vars = "orig.ident")
LAM123_Integrate <- RunUMAP(LAM123_Integrate,  dims = 1:15, reduction = "harmony") |> 
  FindNeighbors(reduction = "harmony", dims = 1:15) |> 
  FindClusters(resolution = 0.5) |> 
  RunTSNE(dims=1:15,reduction = "harmony")

p1 <- DimPlot(LAM123_Integrate, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(LAM123_Integrate, reduction = "tsne", group.by = "seurat_clusters", label = TRUE,repel = TRUE)
p1 + p2
Idents(LAM123_Integrate) = 'seurat_clusters'

##找marker----------------------------------------------------------------------
all.markers <- FindAllMarkers(
  object = LAM123_Integrate,
  only.pos = F, #保留上调与下调的基因
  test.use = 'wilcox', 
  slot = 'data', #默认使用data而不是counts
  min.pct = 0.25, #该基因在至少25%的细胞内表达
  logfc.threshold = 0.25 
) 
#save(all.markers,file = 'ALL_MARKERS.rda')


lam_markers <- c("PMEL", "MLANA", "ACTA2", "DES", "FIGF", "ESR1", "PGR", "CTSK")
present_markers <- lam_markers[lam_markers %in% all.markers$gene]
cat("存在于差异分析结果中的marker基因：", present_markers, "\n")

marker_clusters <- all.markers[all.markers$gene %in% present_markers,]
marker_clusters <- marker_clusters[order(marker_clusters$gene, marker_clusters$cluster), ]
print(marker_clusters)

for (marker in present_markers) {
  p <- VlnPlot(object = LAM123_Integrate, features = marker, group.by = "seurat_clusters") +
    ggtitle(paste("Expression of", marker))
  print(p)
}
DotPlot(object = LAM123_Integrate, features = present_markers) +
  RotatedAxis()

lam_clusters <- 13

# 在元数据中添加细胞类型标签
LAM123_Integrate@meta.data$category <- ifelse(
  LAM123_Integrate@meta.data$seurat_clusters %in% lam_clusters,
  "LAM_cell",
  "Normal_cell"
)

p3 <- DimPlot(LAM123_Integrate, reduction = "tsne", group.by = "category")
p3 + p2

top20 <- all.markers %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>% #以metadata中的‘cluster’来进行分类
  filter(p_val_adj < 0.05) %>%
  top_n(20,avg_log2FC) %>%
  group_by()

library(readxl)
marker <- read_xlsx('Cell_marker_Human.xlsx')
markers <- marker$marker

good_markers <- all.markers %>%
  filter(cluster == 0) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(50) 
cellmarker <- good_markers$gene
cluster0 <- marker[which(marker$marker %in% cellmarker),]
table(cluster0$cell_name) # Alveolar epithelial cell Type 2	
# Astrocyte,Dendritic cell,Epithelial cell

good_markers <- all.markers %>%
  filter(cluster == 1) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(40)
cellmarker <- good_markers$gene
cluster1 <- marker[which(marker$marker %in% cellmarker),]
table(cluster1$cell_name) # Macrophage

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
table(cluster2$cell_name) # Monocyte


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
table(cluster3$cell_name) # Macrophage
#

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
table(cluster4$cell_name) # T cell

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
table(cluster5$cell_name) # Monocyte

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
table(cluster6$cell_name) # Epithelial cell
# Astrocyte 

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
table(cluster7$cell_name) # Dendritic cell ?

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
table(cluster8$cell_name) # ?

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
table(cluster9$cell_name) # Endothelial cell

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
table(cluster10$cell_name) # Ciliated cell

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
table(cluster11$cell_name) # Endothelial cell

good_markers <- all.markers %>%
  filter(cluster == 12) %>%
  filter(avg_log2FC > 0.5,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster12 <- marker[which(marker$marker %in% cellmarker),]
table(cluster12$cell_name) # Alveolar epithelial cell Type 2

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
table(cluster13$cell_name) # LAM cell
#Fibroblast Myofibroblast, Smooth muscle cell

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
table(cluster14$cell_name) # T cell

good_markers <- all.markers %>%
  filter(cluster == 15) %>%
  filter(avg_log2FC > 1,          # 表达量足够高
         p_val_adj < 0.01,        # 统计显著
         pct.1 > 0.5,             # 在该cluster中表达比例高
         pct.2 < 0.3) %>%         # 在其他cluster中表达比例低
  arrange(desc(avg_log2FC)) %>%
  head(20)
cellmarker <- good_markers$gene
cluster15 <- marker[which(marker$marker %in% cellmarker),]
table(cluster15$cell_name) # Plasma cell
# (Memory B cell,Plasmablast)

LAM123_Integrate$cell_type = case_when(
  LAM123_Integrate$seurat_clusters %in% c(2,5) ~ 'Monocyte',
  LAM123_Integrate$seurat_clusters %in% c(4,14) ~ 'T cell',
  LAM123_Integrate$seurat_clusters %in% c(9,11) ~ 'Endothelial cell',
  LAM123_Integrate$seurat_clusters %in% c(6) ~ 'Epithelial cell',
  LAM123_Integrate$seurat_clusters %in% c(10) ~ 'Ciliated cell',
  LAM123_Integrate$seurat_clusters %in% c(8) ~ 'Unknown',
  LAM123_Integrate$seurat_clusters %in% c(1,3) ~ 'Macrophage',
  LAM123_Integrate$seurat_clusters %in% c(13) ~ 'LAM cell',
  LAM123_Integrate$seurat_clusters %in% c(15) ~ 'Plasma cell',
  LAM123_Integrate$seurat_clusters %in% c(7) ~ 'Dendritic cell',
  LAM123_Integrate$seurat_clusters %in% c(0,12) ~ 'Alveolar epithelial cell Type 2',
  TRUE ~ 'Unknown'
)

p4 <- DimPlot(LAM123_Integrate, reduction = "tsne", group.by = "cell_type", label = T, repel = T)
p3 + p4

plotdata = LAM123_Integrate@meta.data
ggplot(
  plotdata,
  aes(x=cell_type,fill = orig.ident)
)+
  geom_bar(position = 'dodge')+  #并排显示分组条形图
  scale_fill_igv(alpha = 0.7) +
  theme(axis.title = element_text(face='bold',hjust = 0.5,size = 16),
        axis.text.y= element_text(face='bold',hjust = 0.5,size = 12)) +
  coord_flip()

# 百分比堆叠图
ggplot(
  plotdata,
  aes(x=orig.ident,fill = cell_type)
)+
  geom_bar(position = 'fill')+ #百分比堆叠
  scale_fill_igv(alpha = 0.7)+
  theme(axis.title = element_text(face='bold',hjust = 0.5,size = 16),
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

Idents(LAM123_Integrate) <- "cell_type"
# 对细胞数多的组进行亚采样
set.seed(123) # 保证结果可重复
macrophage_cells <- WhichCells(LAM123_Integrate, idents = "Macrophage")
subsampled_macrophage <- sample(macrophage_cells, size = 200) # 与LAM cell组数量匹配

# 提取亚采样后的对象进行分析
subset_obj <- subset(LAM123_Integrate, cells = c(WhichCells(LAM123_Integrate, idents = "LAM cell"), 
                                                 subsampled_macrophage))

degs = FindMarkers(
  subset_obj,
  logfc.threshold = 0.5,
  only.pos = F,
  ident.1 = 'LAM cell',
  ident.2 = 'Macrophage',#默认使用ident1 vs ident2，而不是相反
  group.by = 'cell_type'
) %>%
  mutate(gene=rownames(.))

change = as.factor(ifelse(degs$p_val_adj < 0.05 & abs(degs$avg_log2FC) > 1,
                          ifelse(degs$avg_log2FC > 1 ,'Up','Down'),'No change'))#标记上调下调基因
filtered_degs = degs %>%
  filter(p_val_adj<0.01)%>%
  filter(abs(avg_log2FC)>5)
filtered_degs <- filtered_degs[1:40,]
degs$label <- ifelse(degs$gene %in% filtered_degs$gene, degs$gene, NA)#标记相当显著的基因

ggplot(degs, aes(avg_log2FC, -log10(p_val_adj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = label), 
                  size = 3, 
                  box.padding = 0.2, 
                  max.overlaps = 10) + # 防止重叠
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No change" = "gray"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right')+
  scale_x_continuous(limits = c(-max(abs(degs$avg_log2FC)), max(abs(degs$avg_log2FC))))


##富集分析----------------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(enrichplot)

LAM_degs = all.markers %>%
  filter(cluster == 13) %>%
  filter(p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 0.25) %>%
  arrange(desc(avg_log2FC))

GO_database <- 'org.Hs.eg.db'
ENTgene <- bitr(LAM_degs$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database,drop = F) %>%
  distinct(SYMBOL,.keep_all = T) %>%
  drop_na()
ENTgene <- dplyr::distinct(ENTgene,SYMBOL,.keep_all = T)
GO <- enrichGO(
  ENTgene$ENTREZID,
  OrgDb = GO_database, #即'org.Hs.eg.db'，即人类数据库
  keyType = "ENTREZID",
  ont = "all",  # 获取BP, CC, MF所有本体论的富集结果
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.01,
  readable = T
)
barplot(
  GO,
  split="ONTOLOGY", #按基因本体论分组
)+
  facet_wrap(~ONTOLOGY, scales="free", ncol=1)
dotplot(GO,split = "ONTOLOGY")

KEGG_database = 'hsa'
kegg_data <- enrichKEGG(
  gene = ENTgene$ENTREZID,
  keyType = 'kegg',
  organism = KEGG_database,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  use_internal_data = F
) %>%
  setReadable(.,OrgDb = GO_database,keyType = 'ENTREZID')

kegg_data = setReadable(kegg_data, #前面分析的结果
                        OrgDb = "org.Hs.eg.db", #人类数据库
                        keyType = "ENTREZID") 
barplot(kegg_data, x = "GeneRatio", color = "p.adjust", #默认参数
        showCategory =15)
dotplot(kegg_data,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory = 15) 

LAM_degs$SYMBOL = LAM_degs$gene
cnet_data = merge(LAM_degs,ENTgene,by = "SYMBOL")
cnet_data = data.frame(cnet_data$ENTREZID,cnet_data$gene,cnet_data$avg_log2FC)
logFC <- cnet_data$cnet_data.avg_log2FC
names(logFC) <- cnet_data$cnet_data.gene
cnetplot(
  kegg_data, 
  showCategory = 8,  # 选择top8的pathway
  foldChange = logFC,
  colorEdge = TRUE,
  node_label = 'all',
  max.overlaps = 50
) + 
  # 三色渐变：低表达(蓝色) -> 中间值(灰色) -> 高表达(红色)
  scale_color_gradient2(
    low = "#3366CC",    # 低fold change颜色（蓝色）
    mid = "gray",       # 中间值颜色（灰色）
    high = "#DC3912",   # 高fold change颜色（红色）
    midpoint = 0        # 中间值点（fold change=0处）
  )

##看通路--------------------------------------------------------------------------
library(pathview)
pathview_data = LAM_degs %>%
  filter(LAM_degs$gene %in% as.factor(ENTgene$SYMBOL))
pathview_data$SYMBOL = pathview_data$gene
pathview_data = merge(pathview_data,ENTgene,by = 'SYMBOL')
pathview_data_draw = data.frame(pathview_data$ENTREZID,pathview_data$avg_log2FC) 
rownames(pathview_data_draw) = pathview_data_draw$pathview_data.ENTREZID
colnames(pathview_data_draw) <- c('ENTREZID','log2FC')
pathview_data_draw$ENTREZID <- as.numeric(pathview_data_draw$ENTREZID)

# 转换为pathview需要的格式：命名数值向量（名称为ENTREZID，值为log2FC）
pathview_vector <- pathview_data$avg_log2FC
names(pathview_vector) <- pathview_data$ENTREZID

# 拆分上调和下调基因
up_genes <- pathview_vector[pathview_vector > 0]
down_genes <- pathview_vector[pathview_vector < 0] 

pathview(gene.data = up_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa05152", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Tuberculosis_up", #输出文件名
         kegg.native = T
)
pathview(gene.data = down_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa05152", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Tuberculosis_down", #输出文件名
         kegg.native = T
)

#save(LAM123_Integrate, file = 'LAM123_Integrate.rda')

