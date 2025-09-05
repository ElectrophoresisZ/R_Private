files <- dir(pattern = 'gz$')
sapply(files, R.utils::gunzip)
files <- dir(pattern = 'csv$')
for (i in 1:length(files)){
  file.rename(files[i], paste(strsplit(files[i],'_')[[1]][1],'.csv',sep = ''))
}
files <- dir(pattern = 'csv$')

library(tidyverse)
library(readr)

count_matrix <- data.frame(Response=rep(0,18446))
for (i in 1:length(files)){
  data <- read.csv(files[i], sep = '\t', header = T)
  #Gene_Symbol <- gsub('\\.csv$','',files[i])
  count_matrix <- cbind(count_matrix, data$tpm)
}

count_matrix <- do.call(cbind, lapply(files, function(file){
  data <- read.csv(file, sep = '\t')[, c('target_id','tpm')]
  colnames(data) <- c('Genes', gsub('\\.csv$','',file))
  return(data)
}))
ID_set <- gsub('\\.csv$','',files)

rawdata <- count_matrix[, grep('tpm',colnames(count_matrix))]
colnames(rawdata) <- ID_set
gene_set <- rownames(rawdata)
for (i in 1:length(gene_set)){
  gene_set[i] <- as.character(strsplit(gene_set[i],'_')[[1]][1])
}
rownames(rawdata) <- gene_set 

library(edgeR)
dge <- DGEList(rawdata)
dge <- calcNormFactors(dge, method = 'TMM')
exp_GSE152728 <- dge$counts
#save(exp_GSE152728, file = 'exp_GSE152728.Rdata')


library(GEOquery)
gset_GSE152728 <- getGEO('GSE152728', getGPL = F)
Gset_GSE152728 <- gset_GSE152728[[1]]
pdata_GSE152728 <- pData(Gset_GSE152728)
table(pdata_GSE152728$characteristics_ch1.1)
#save(pdata_GSE152728,file = 'pdata_GSE152728.Rdata')
