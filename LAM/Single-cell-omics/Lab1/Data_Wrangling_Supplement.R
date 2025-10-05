library(tidyverse)
library(clusterProfiler)
library(enrichplot)

load('LAM_Macrophage_degs_enrich.rda')
LAM_Macrophage_degs_enrich_up <- LAM_Macrophage_degs_enrich %>%
  filter(avg_log2FC>0)

LAM_Macrophage_degs_enrich_down <- LAM_Macrophage_degs_enrich %>%
  filter(avg_log2FC<0)

##上调--------------------------------------------------------------------------
GO_database <- 'org.Hs.eg.db'
ENTgene <- bitr(LAM_Macrophage_degs_enrich_up$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database,drop = F) %>%
  distinct(SYMBOL,.keep_all = T) %>%
  drop_na()
ENTgene <- dplyr::distinct(ENTgene,SYMBOL,.keep_all = T)
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

LAM_Macrophage_degs_enrich_up$SYMBOL = LAM_Macrophage_degs_enrich_up$gene
cnet_data = merge(LAM_Macrophage_degs_enrich_up,ENTgene,by = "SYMBOL")
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

library(pathview)
pathview_data = LAM_Macrophage_degs_enrich_up %>%
  filter(LAM_Macrophage_degs_enrich_up$gene %in% as.factor(ENTgene$SYMBOL))
pathview_data$SYMBOL = pathview_data$gene
pathview_data = merge(pathview_data,ENTgene,by = 'SYMBOL')
pathview_data_draw = data.frame(pathview_data$ENTREZID,pathview_data$avg_log2FC) 
rownames(pathview_data_draw) = pathview_data_draw$pathview_data.ENTREZID
colnames(pathview_data_draw) <- c('ENTREZID','log2FC')
pathview_data_draw$ENTREZID <- as.numeric(pathview_data_draw$ENTREZID)

# 转换为pathview需要的格式：命名数值向量（名称为ENTREZID，值为log2FC）
pathview_vector <- pathview_data$avg_log2FC
names(pathview_vector) <- pathview_data$ENTREZID
#save(pathview_data, file = 'LAM_Macrophage_pathview_up.rda')

# 拆分上调和下调基因
up_genes <- pathview_vector[pathview_vector > 0]
down_genes <- pathview_vector[pathview_vector < 0] 

pathview(gene.data = up_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04151", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "PI3K-AKT_Signaling_Pathway_up", #输出文件名
         kegg.native = T
)
pathview(gene.data = down_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04151", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "PI3K-AKT_Signaling_Pathway_down", #输出文件名
         kegg.native = T
)

pathview(gene.data = up_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04810", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Regulation_of_Actin_Cytoskeleton_up", #输出文件名
         kegg.native = T
)
pathview(gene.data = down_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04810", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Regulation_of_Actin_Cytoskeleton_down", #输出文件名
         kegg.native = T
)


##下调--------------------------------------------------------------------------
ENTgene <- bitr(LAM_Macrophage_degs_enrich_down$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database,drop = F) %>%
  distinct(SYMBOL,.keep_all = T) %>%
  drop_na()
ENTgene <- dplyr::distinct(ENTgene,SYMBOL,.keep_all = T)
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

LAM_Macrophage_degs_enrich_down$SYMBOL = LAM_Macrophage_degs_enrich_down$gene
cnet_data = merge(LAM_Macrophage_degs_enrich_down,ENTgene,by = "SYMBOL")
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

library(pathview)
pathview_data = LAM_Macrophage_degs_enrich_down %>%
  filter(LAM_Macrophage_degs_enrich_down$gene %in% as.factor(ENTgene$SYMBOL))
pathview_data$SYMBOL = pathview_data$gene
pathview_data = merge(pathview_data,ENTgene,by = 'SYMBOL')
pathview_data_draw = data.frame(pathview_data$ENTREZID,pathview_data$avg_log2FC) 
rownames(pathview_data_draw) = pathview_data_draw$pathview_data.ENTREZID
colnames(pathview_data_draw) <- c('ENTREZID','log2FC')
pathview_data_draw$ENTREZID <- as.numeric(pathview_data_draw$ENTREZID)

# 转换为pathview需要的格式：命名数值向量（名称为ENTREZID，值为log2FC）
pathview_vector <- pathview_data$avg_log2FC
names(pathview_vector) <- pathview_data$ENTREZID
#save(pathview_data, file = 'LAM_Macrophage_pathview_down.rda')

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

pathview(gene.data = up_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04612", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Antigen_Processing_and_Presentation_up", #输出文件名
         kegg.native = T
)
pathview(gene.data = down_genes, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04612", #选择一个KEGG信号通路,这里选的是hsa04650。更多通路见KEGG官网
         species = "hsa",
         out.suffix = "Antigen_Processing_and_Presentation_down", #输出文件名
         kegg.native = T
)


