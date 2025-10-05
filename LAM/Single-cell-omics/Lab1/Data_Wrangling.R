library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)

##数据读取----------------------------------------------------------------------
LAM_lab1 <- readRDS('seurat.combined.sc.rds')
table(LAM_lab1$new.ident)
table(LAM_lab1$treatment)

p1 <- DimPlot(LAM_lab1, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(LAM_lab1, reduction = "umap", group.by = "seurat.clusters", label = TRUE,repel = TRUE)
p1
p2

##UMAP图对比--------------------------------------------------------------------
LAM_lab1@meta.data$category <- ifelse(
  LAM_lab1@meta.data$new.ident == 'LAM',
  "LAM_cell",
  "Normal_cell"
)

p1 <- DimPlot(LAM_lab1, reduction = "umap", group.by = "treatment")
p3 <- DimPlot(LAM_lab1, reduction = "umap", group.by = "new.ident")
p4 <- DimPlot(LAM_lab1, reduction = "umap", group.by = "category")
p1+p4+p3

##百分比堆叠图------------------------------------------------------------------
plotdata = LAM_lab1@meta.data
ggplot(
  plotdata,
  aes(x=new.ident,fill = treatment)
)+
  geom_bar(position = 'dodge')+  #并排显示分组条形图
  scale_fill_igv(alpha = 0.7) +
  theme(axis.title = element_text(face='bold',hjust = 0.5,size = 16),
        axis.text.y= element_text(face='bold',hjust = 0.5,size = 12)) +
  coord_flip()

ggplot(
  plotdata,
  aes(x=treatment,fill = new.ident)
)+
  geom_bar(position = 'fill')+ #百分比堆叠
  scale_fill_igv(alpha = 0.7)+
  theme(axis.title = element_text(face='bold',hjust = 0.5,size = 16),
        axis.text.x= element_text(face='bold',hjust = 0.5,size = 12))

Idents(LAM_lab1) <- "new.ident"


##对细胞数多的组进行亚采样------------------------------------------------------
table(LAM_lab1$new.ident)
table(LAM_lab1$treatment)

set.seed(123) # 保证结果可重复
LAM_lab1_naive <- subset(LAM_lab1, subset = treatment == 'Naïve')
table(LAM_lab1_naive$treatment)
macrophage_cells <- WhichCells(LAM_lab1_naive, idents = "Macrophages")
subsampled_macrophage <- sample(macrophage_cells, size = 110) # 与LAM cell组数量匹配

# 提取亚采样后的对象进行分析
subset_LAM_Macrophage <- subset(LAM_lab1, cells = c(WhichCells(LAM_lab1, idents = "LAM"), 
                                                 subsampled_macrophage))
table(subset_LAM_Macrophage$treatment)

LAM_Macrophage_degs = FindMarkers(
  subset_LAM_Macrophage,
  logfc.threshold = 0.5,
  only.pos = F,
  ident.1 = 'LAM',
  ident.2 = 'Macrophages',#默认使用ident1 vs ident2，而不是相反
  group.by = 'new.ident'
) %>%
  mutate(gene=rownames(.))

change = as.factor(ifelse(LAM_Macrophage_degs$p_val_adj < 0.05 & abs(LAM_Macrophage_degs$avg_log2FC) > 1,
                          ifelse(LAM_Macrophage_degs$avg_log2FC > 1 ,'Up','Down'),'No change'))#标记上调下调基因
filtered_LAM_Macrophage_degs = LAM_Macrophage_degs %>%
  filter(p_val_adj<0.01)%>%
  filter(abs(avg_log2FC)>8)
filtered_LAM_Macrophage_degs <- filtered_LAM_Macrophage_degs[1:40,]
LAM_Macrophage_degs$label <- ifelse(LAM_Macrophage_degs$gene %in% filtered_LAM_Macrophage_degs$gene, LAM_Macrophage_degs$gene, NA)#标记相当显著的基因

ggplot(LAM_Macrophage_degs, aes(avg_log2FC, -log10(p_val_adj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.5, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = label), 
                  size = 3, 
                  box.padding = 0.2, 
                  max.overlaps = 10) + # 防止重叠
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No change" = "gray"))+
  theme(panel.grid = element_blank(),
        legend.position = 'right')+
  scale_x_continuous(limits = c(-max(abs(LAM_Macrophage_degs$avg_log2FC)), max(abs(LAM_Macrophage_degs$avg_log2FC))))


##富集分析----------------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(enrichplot)

LAM_Macrophage_degs_enrich <- LAM_Macrophage_degs %>%
  filter(p_val_adj<0.05) %>%
  filter(abs(avg_log2FC)>0.25) %>%
  arrange(desc(avg_log2FC))
#save(LAM_Macrophage_degs_enrich,file = 'LAM_Macrophage_degs_enrich.rda')

GO_database <- 'org.Hs.eg.db'
ENTgene <- bitr(LAM_Macrophage_degs_enrich$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database,drop = F) %>%
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

LAM_Macrophage_degs$SYMBOL = LAM_Macrophage_degs$gene
cnet_data = merge(LAM_Macrophage_degs,ENTgene,by = "SYMBOL")
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
pathview_data = LAM_Macrophage_degs %>%
  filter(LAM_Macrophage_degs$gene %in% as.factor(ENTgene$SYMBOL))
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

#save(LAM_lab1, file = 'LAM_lab1.rda')
