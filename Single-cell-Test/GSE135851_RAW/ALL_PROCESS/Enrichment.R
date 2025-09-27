library(tidyverse)
library(clusterProfiler)
library(enrichplot)

load("MANUAL_ALL_MARKERS.rda")
Epithelial_degs = cell.markers %>%
  filter(cluster == "Epithelial cell") %>%
  filter(p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC))

#GO------------------------------------------------------------------------------
GO_database <- 'org.Hs.eg.db'
ENTgene <- bitr(Epithelial_degs$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database,drop = F) %>%   # drop:保留无法匹配的基因
  distinct(SYMBOL,.keep_all = T) %>% #基于基因符号列去重,保留所有其他列的数据
  drop_na() #保留所有其他列的数据
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

#Biological Process (BP)
#Cellular Component (CC)
#Molecular Function (MF)
barplot(
  GO,
  split="ONTOLOGY", #按基因本体论分组
)+
  facet_wrap(~ONTOLOGY, scales="free", ncol=1) #每个面板根据自身数据调整坐标轴
dotplot(GO,split = "ONTOLOGY")

#KEGG------------------------------------------------------------------
KEGG_database = 'hsa'
kegg_data <- enrichKEGG(
  gene = ENTgene$ENTREZID,
  keyType = 'kegg',
  organism = KEGG_database,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  use_internal_data = F #使用在线或本地数据
) %>%
  setReadable(.,OrgDb = GO_database,keyType = 'ENTREZID') #管道输入的数据,将Entrez ID转换为基因符号
kegg_data = setReadable(kegg_data, #前面分析的结果
                        OrgDb = "org.Hs.eg.db", #人类数据库
                        keyType = "ENTREZID") #要转换的基因类型
#save(kegg_data,file = 'kegg_epithelial_data.rda')
#write.table(kegg_data,file="kegg_epithelial_data.txt",sep="\t",quote=F,row.names = F)

barplot(kegg_data, x = "GeneRatio", color = "p.adjust", #默认参数
        showCategory =10) #只显示前10
dotplot(kegg_data,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory = 10) 

Epithelial_degs$SYMBOL = Epithelial_degs$gene
cnet_data = merge(Epithelial_degs,ENTgene,by = "SYMBOL")
cnet_data = data.frame(cnet_data$ENTREZID,cnet_data$gene,cnet_data$avg_log2FC)
logFC <- cnet_data$cnet_data.avg_log2FC
names(logFC) <- cnet_data$cnet_data.gene

cnetplot(kegg_data, showCategory = 8, #选择top8的pathway ，这里也可以用包含pathway名称的向量           
         foldChange = logFC, #基因表达变化数据
         colorEdge = T,
         node_label = 'all', #显示所有节点的标签（基因和通路名称
         max.overlaps = 50 #最多允许50个标签重叠，超过时自动隐藏部分标签
)

#GSEA-------------------------------------------------------------------------------
expr_matrix = data.frame(Epithelial_degs$SYMBOL,Epithelial_degs$avg_log2FC)
colnames(expr_matrix) <- c("SYMBOL","Log2FoldChange")
merged_matrix = merge(expr_matrix,ENTgene,by = "SYMBOL")
data_sort <- arrange(merged_matrix,desc(Log2FoldChange)) #log2FC进行排序，这是GSEA分析必要的
gene_list <- data_sort$Log2FoldChange
names(gene_list) <- data_sort$ENTREZID

KEGG_database = 'hsa'
gsea <- gseKEGG(gene_list,organism = KEGG_database, pvalueCutoff = 0.05)
gsea = setReadable(gsea, OrgDb = GO_database, keyType = "ENTREZID")

dotplot(gsea)
ridgeplot(gsea,label_format = 100)
gseaplot2(gsea,c(1:10),pvalue_table = F)

#---------------------------------------------------------------------------------
library(pathview)
library(tidyverse)

#在这张现成的、标准化的通路上，用颜色标记出哪些（基因/蛋白）目前正处于上调或下调状态
pathview_data = Epithelial_degs %>%
  filter(Epithelial_degs$gene %in% as.factor(ENTgene$SYMBOL))
pathview_data$SYMBOL = pathview_data$gene
pathview_data = merge(pathview_data,ENTgene,by = 'SYMBOL')
pathview_data_draw = data.frame(pathview_data$ENTREZID,pathview_data$avg_log2FC) 
rownames(pathview_data_draw) = pathview_data_draw$pathview_data.ENTREZID
colnames(pathview_data_draw) <- c('ENTREZID','log2FC')
pathview_data_draw$ENTREZID <- as.numeric(pathview_data_draw$ENTREZID)

pathview(gene.data = pathview_data_draw, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "05152", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Tuberculosis_Pathway", #输出文件名
         kegg.native = T
)
pathview(gene.data = pathview_data_draw, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "04145", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Phagosome_Pathway", #输出文件名
         kegg.native = T
)
