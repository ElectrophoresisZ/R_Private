load('macrophage_gene.Rdata')
library(clusterProfiler)
library(org.Hs.eg.db)
mRNA <- bitr(rownames(macrophage_gene),
             fromType = 'SYMBOL', toType = 'ENTREZID',
             OrgDb = 'org.Hs.eg.db')
rownames(mRNA) <- mRNA$ENTREZID
gene_trans <- macrophage_gene[mRNA$SYMBOL,]
macro_gsea <- as.data.frame(cbind(gene_trans[,2],mRNA[,2]))
colnames(macro_gsea) <- c('avglog2FC','ENTREZID')
rownames(macro_gsea) <- mRNA$SYMBOL
macro_gsea <- macro_gsea[order(as.numeric(macro_gsea$avglog2FC),decreasing = T),]
head(macro_gsea$avglog2FC)
#save(macro_gsea,file = 'Macrophage_GSEA.Rdata')
macro_gsea_fc <- as.numeric(macro_gsea$avglog2FC)
names(macro_gsea_fc) <- macro_gsea$ENTREZID

macro_gsea_result <- gseKEGG(macro_gsea_fc,organism = 'hsa',keyType = 'kegg',
                             minGSSize = 10, maxGSSize = 500,pvalueCutoff = 0.05,
                             pAdjustMethod = 'BH')
dim(macro_gsea_result)
#save(macro_gsea_result,file = 'Macro_GSEA_result.Rdata')

library(DOSE)
dotplot(macro_gsea_result, showCategory=15)


library(enrichplot)







