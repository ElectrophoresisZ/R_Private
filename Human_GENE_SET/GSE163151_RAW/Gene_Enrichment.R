rm(list=ls())
load('GSE163151_updated.Rdata')
load('perm.pVal.Rdata')
selected_gene <- read.table('GSE163151_select.txt')


GSE163151_genes <- selected_gene

library(clusterProfiler)
library(org.Hs.eg.db)
GSE163151_genes_mRNA <- bitr(row.names(GSE163151_genes), 
                             fromType = 'SYMBOL', toType = 'ENTREZID',
                             OrgDb = 'org.Hs.eg.db')
rownames(GSE163151_genes_mRNA) <- GSE163151_genes_mRNA$ENTREZID
GSE163151_genes_analysis <- GSE163151_genes[GSE163151_genes_mRNA[,1],]
GSE163151_GSEA <- as.data.frame(cbind(GSE163151_genes_mRNA[,2], 
                                      as.numeric(GSE163151_genes_analysis[,1])))
colnames(GSE163151_GSEA) <- c('ENTREZID','logFC')
GSE163151_GSEA <- GSE163151_GSEA[order(GSE163151_GSEA$logFC, decreasing = T),]
row.names(GSE163151_GSEA) <- GSE163151_GSEA[,1]
GSE163151_GSEA_FC <- sort(as.numeric(GSE163151_GSEA[,2]), decreasing = T)
names(GSE163151_GSEA_FC) <- rownames(GSE163151_GSEA)
save(GSE163151_GSEA, file = 'GSEA_Mapping.Rdata')

gsea_GSE163151 <- gseKEGG(GSE163151_GSEA_FC, organism = 'hsa',
                          keyType = 'kegg', minGSSize = 10, maxGSSize = 500,
                          pvalueCutoff = 0.05, pAdjustMethod = 'BH')
save(gsea_GSE163151, file = 'GSEA_Analysis_GSE163151.Rdata')
dim(gsea_GSE163151)

library(DOSE)
dotplot(gsea_GSE163151, showCategory=15)

library(enrichplot)

gseaplot2(gsea_GSE163151, row.names(data.frame(gsea_GSE163151))[1],
          title = data.frame(gsea_GSE163151)[1,'Description'], pvalue_table = F)
gseaplot2(gsea_GSE163151, row.names(data.frame(gsea_GSE163151))[1:5], pvalue_table = F)

# q值 (FDR)	假发现率（False Discovery Rate, FDR）
gsea_core <- data.frame(gsea_GSE163151)
gsea_core$genes <- gsea_core$core_enrichment
dim(gsea_core)
gsea_genelist <- strsplit(gsea_core$core_enrichment,'/')
genes_SYMBOL <- c()

for (i in 1:nrow(gsea_core)){
  genes_ENTREZID <- gsea_genelist[[i]]
  for (j in genes_ENTREZID){
    genes_SYMBOL <- c(genes_SYMBOL, GSE163151_genes_mRNA[j,1])
  }
  gsea_core[i,12] <- paste(genes_SYMBOL, collapse = '/')
  genes_SYMBOL <- c()
}

# 不会把数据框的行名写进文件里
# 不会给字符型数据加双引号
write.table(gsea_core, file = 'GSEA_table.txt', sep = '\t', row.names = F, quote=F)

#-----------------------------------------------------------------------------

# 另一种富集分析方法
library(fgsea)
library(doParallel)
library(RCPA)
kegg_geneSets <- getGeneSets(org = "hsa", database = "KEGG")
fgsea_GSE163151 <- fgseaMultilevel(pathways = kegg_geneSets,
                                   stats = GSE163151_GSEA_FC,
                                   minSize = 15,
                                   maxSize = 500,
                                   eps = 0.0,
                                   nproc = 4)






