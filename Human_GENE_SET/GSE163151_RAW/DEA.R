rm(list = ls())
load('exp_data_subset_GSE163151.Rdata')
exp_GSE163151 <- log2_cpm_subset
load('pdata_GSE163151_subset.Rdata')
table(pdata_GSE163151_subset$`disease state:ch1`)
# COVID-19                 Non-viral acute respiratory illness 
#    145                                  82
sign <- pdata_GSE163151_subset$`disease state:ch1`
design <- model.matrix(~0+factor(sign)) # 不要截距项
colnames(design) <- c('COVID.19','Non.ARI')

library(limma)
fit <- lmFit(exp_GSE163151,design)
cont.matrix <- makeContrasts(COVID.19-Non.ARI, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEA_GSE163151 <- topTable(fit2,coef = 1, number = dim(exp_GSE163151)[1],adjust.method = 'BH')
#write.table(DEA_GSE163151, file = 'DEA_GSE163151.txt',sep = '\t',quote = F)

GSE163151_selectP <- DEA_GSE163151[which(as.numeric(DEA_GSE163151$adj.P.Val)<0.05),]
dim(GSE163151_selectP)
GSE163151_selectUp <- GSE163151_selectP[which(as.numeric(GSE163151_selectP$logFC)>log2(3/2)),]
dim(GSE163151_selectUp)
GSE163151_selectDown <- GSE163151_selectP[which(as.numeric(GSE163151_selectP$logFC)<log2(2/3)),]
dim(GSE163151_selectDown)

GSE163151_select_final <- rbind(GSE163151_selectUp,GSE163151_selectDown)
GSE163151_select_final <- GSE163151_select_final[order(as.numeric(GSE163151_select_final$logFC), decreasing = T),]
#write.table(GSE163151_select_final, file='GSE163151_select.txt',sep = '\t',quote = F)

library(ggplot2)
DEA_GSE163151$sign <- 'no'
DEA_GSE163151[row.names(GSE163151_selectUp),7] <- 'up'
DEA_GSE163151[row.names(GSE163151_selectDown),7] <- 'down'
ggplot(DEA_GSE163151, aes(x=logFC,y=-log10(as.numeric(P.Value)),colour = sign)) +
  geom_point() +
  scale_color_manual(values=c("royalblue4","grey70","firebrick3"))

library(pheatmap)  
name_select <- row.names(rbind(GSE163151_select_final[1:20,], GSE163151_select_final[7084:7103,]))
GSE163151_plot <- exp_GSE163151[name_select,]
annotation_col <- data.frame(Type = factor(sign))
row.names(annotation_col) <- colnames(exp_GSE163151)
pheatmap(GSE163151_plot, scale = 'row',
         show_rownames = T, show_colnames = F,
         cluster_rows = F, cluster_cols = T,
         annotation = annotation_col,
         color = colorRampPalette(c("royalblue4","white","firebrick3"))(500))

#save(GSE163151_plot, file = '40selected_genes.Rdata')

name_select2 <- row.names(rbind(GSE163151_selectUp[1:100,],GSE163151_selectDown[1:100,]))
GSE163151_genes <- exp_GSE163151[name_select2,]
#save(GSE163151_genes,file = '200selected_genes.Rdata')
