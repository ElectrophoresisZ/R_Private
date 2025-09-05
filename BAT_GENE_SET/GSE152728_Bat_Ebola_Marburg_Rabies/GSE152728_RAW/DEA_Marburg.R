rm(list = ls())
load('exp_GSE152728.Rdata')
load('pdata_GSE152728.Rdata')

exp_GSE152728 <- as.data.frame(exp_GSE152728)
table(pdata_GSE152728$characteristics_ch1.1)

pdata_subset_Marburg <-pdata_GSE152728[which(pdata_GSE152728$characteristics_ch1.1 %in%
                                             c('treatment: infected with Marburg','treatment: infected with NAIVE')),]
#save(pdata_subset_Marburg, file = 'pdata_subset_Marburg.Rdata')
exp_subset_Marburg <- exp_GSE152728[,rownames(pdata_subset_Marburg)]

sign <- pdata_subset_Marburg$characteristics_ch1.1
design_matrix <- model.matrix(~0+factor(sign, 
                                        levels=c('treatment: infected with NAIVE','treatment: infected with Marburg')))
colnames(design_matrix) <- c('NAIVE','Marburg')

library(limma)
v_GSE152728 <- voom(exp_subset_Marburg, design_matrix)
exp_subset_Marburg <- v_GSE152728$E
#save(exp_subset_Marburg,file = 'exp_subset_Marburg.Rdata')
cont.matrix <- makeContrasts(Marburg-NAIVE, levels = design_matrix)
fit <- lmFit(exp_subset_Marburg, design_matrix)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEA_Marburg <- topTable(fit2,coef=1,number=dim(exp_subset_Marburg)[1],adjust.method = 'BH')

DEA_Marburg_SelectP <- DEA_Marburg[which(as.numeric(DEA_Marburg$adj.P.Val)<0.05),]
DEA_Marburg_Up <- DEA_Marburg_SelectP[which(as.numeric(DEA_Marburg_SelectP$logFC)>log2(3/2)),]
dim(DEA_Marburg_Up)
DEA_Marburg_Down <- DEA_Marburg_SelectP[which(as.numeric(DEA_Marburg_SelectP$logFC)<log2(2/3)),]
dim(DEA_Marburg_Down)

DEA_Marburg_Final <- rbind(DEA_Marburg_Up,DEA_Marburg_Down)
DEA_Marburg_Final <- DEA_Marburg_Final[order(as.numeric(DEA_Marburg_Final$logFC),decreasing = T),]
#write.table(DEA_Marburg_Final, file = 'DEA_Marburg_Final.txt', sep = '\t', quote = F)

library(pheatmap)
nrow(DEA_Marburg_Final)
name_select <- rownames(rbind(DEA_Marburg_Final[1:20,],DEA_Marburg_Final[4636:4655,]))
Marburg_plot <- exp_subset_Marburg[name_select,]
annotation_col <- data.frame(Type=factor(sign, 
                                         levels=c('treatment: infected with NAIVE','treatment: infected with Marburg')))
rownames(annotation_col) <- colnames(Marburg_plot)
pheatmap(Marburg_plot, scale = 'row',
         show_rownames = T, show_colnames = F,
         cluster_rows = F, cluster_cols = T,
         annotation_col = annotation_col,
         color = colorRampPalette(c('royalblue','white','red'))(500))
