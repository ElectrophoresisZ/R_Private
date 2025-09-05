rm(list = ls())
load('exp_GSE152728.Rdata')
load('pdata_GSE152728.Rdata')

exp_GSE152728 <- as.data.frame(exp_GSE152728)
table(pdata_GSE152728$characteristics_ch1.1)

pdata_subset_Ebola <-pdata_GSE152728[which(pdata_GSE152728$characteristics_ch1.1 %in%
                                             c('treatment: infected with Ebola','treatment: infected with NAIVE')),]
#save(pdata_subset_Ebola, file = 'pdata_subset_Ebola.Rdata')
exp_subset_Ebola <- exp_GSE152728[,rownames(pdata_subset_Ebola)]

sign <- pdata_subset_Ebola$characteristics_ch1.1
design_matrix <- model.matrix(~0+factor(sign, 
                                        levels=c('treatment: infected with NAIVE','treatment: infected with Ebola')))
colnames(design_matrix) <- c('NAIVE','Ebola')

library(limma)
v_GSE152728 <- voom(exp_subset_Ebola, design_matrix)
exp_subset_Ebola <- v_GSE152728$E
#save(exp_subset_Ebola,file = 'exp_subset_Ebola.Rdata')
cont.matrix <- makeContrasts(Ebola-NAIVE, levels = design_matrix)
fit <- lmFit(exp_subset_Ebola, design_matrix)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEA_Ebola <- topTable(fit2,coef=1,number=dim(exp_subset_Ebola)[1],adjust.method = 'BH')

DEA_Ebola_SelectP <- DEA_Ebola[which(as.numeric(DEA_Ebola$adj.P.Val)<0.05),]
DEA_Ebola_Up <- DEA_Ebola_SelectP[which(as.numeric(DEA_Ebola_SelectP$logFC)>log2(3/2)),]
dim(DEA_Ebola_Up)
DEA_Ebola_Down <- DEA_Ebola_SelectP[which(as.numeric(DEA_Ebola_SelectP$logFC)<log2(2/3)),]
dim(DEA_Ebola_Down)

DEA_Ebola_Final <- rbind(DEA_Ebola_Up,DEA_Ebola_Down)
DEA_Ebola_Final <- DEA_Ebola_Final[order(as.numeric(DEA_Ebola_Final$logFC),decreasing = T),]
#write.table(DEA_Ebola_Final, file = 'DEA_Ebola_Final.txt', sep = '\t', quote = F)

library(pheatmap)
nrow(DEA_Ebola_Final)
name_select <- rownames(rbind(DEA_Ebola_Final[1:20,],DEA_Ebola_Final[1765:1784,]))
Ebola_plot <- exp_subset_Ebola[name_select,]
annotation_col <- data.frame(Type=factor(sign, 
                                         levels=c('treatment: infected with NAIVE','treatment: infected with Ebola')))
rownames(annotation_col) <- colnames(Ebola_plot)
pheatmap(Ebola_plot, scale = 'row',
         show_rownames = T, show_colnames = F,
         cluster_rows = F, cluster_cols = T,
         annotation_col = annotation_col,
         color = colorRampPalette(c('royalblue','white','red'))(500))
