files <- dir(pattern = 'gz$')
sapply(files, R.utils::gunzip)
cellfiles <- dir(pattern = 'txt$')
for (i in 1:length(cellfiles)){
  file.rename(cellfiles[i], paste(strsplit(cellfiles[i],'_')[[1]][1], '.txt', sep=''))
}
cellfiles <- dir(pattern = 'txt$')

library(tidyverse)
count_matrix <- cellfiles %>%
  set_names(gsub('\\.txt$','',cellfiles)) %>%
  map_dfc(~ read.table(.x, header = T, row.names = 1)) # 将文件的第一列作为行名
head(count_matrix)

library(edgeR)
dge <- DGEList(counts = count_matrix)

dge <- calcNormFactors(dge, method = 'TMM')
log2_cpm <- cpm(dge, log = T, prior.count = 1)
colnames(log2_cpm) <- unlist(lapply(cellfiles,strsplit,'\\.'))[c(1:length(cellfiles)*2-1)]
save(log2_cpm, file = 'exp_data_GSE163151.Rdata')

hist(log2_cpm, main = 'rma')
par(mfrow = c(1,1),cex.lab=0.7,cex.axis=0.7)
library("RColorBrewer") 
cols <- brewer.pal(8, "Set1")
boxplot(log2_cpm, lwd=1, col=cols, las=2, 
        names=unlist(lapply(cellfiles,strsplit,'\\.'))[c(1:length(cellfiles)*2-1)])

library(GEOquery)
gset_GSE163151 <- getGEO('GSE163151',getGPL = F)
save(gset_GSE163151,file = 'gset_GSE163151.Rdata')

Gset_GSE163151 <- gset_GSE163151[[1]]
pdata_GSE163151 <- pData(Gset_GSE163151)
table(pdata_GSE163151$`disease state:ch1`)
pdata_GSE163151_subset <- pdata_GSE163151[pdata_GSE163151$`disease state:ch1` %in% c('COVID-19','Non-viral acute respiratory illness'),]
save(pdata_GSE163151, file = 'pdata_GSE163151.Rdata')
save(pdata_GSE163151_subset, file= 'pdata_GSE163151_subset.Rdata')

sample_subset <- as.vector(row.names(pdata_GSE163151_subset))
log2_cpm_subset <- log2_cpm[,colnames(log2_cpm) %in% sample_subset]
save(log2_cpm_subset,file = 'exp_data_subset_GSE163151.Rdata')

plat_GPL24676 <- getGEO('GPL24676')
save(plat_GPL24676, file = 'plat_GPL24676.Rdata')

du <- unique(log2_cpm_subset[duplicated(rownames(log2_cpm_subset))])
length(du)

