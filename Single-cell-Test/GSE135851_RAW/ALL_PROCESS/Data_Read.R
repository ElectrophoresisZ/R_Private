library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)

LAM1.data <- Read10X('LAM1/')
hist(log2(LAM1.data@x),breaks=100)
t_gene <- dimnames(LAM1.data)[[1]]
t_cell <- dimnames(LAM1.data)[[2]]
metadata_t <- data.frame(row.names = t_cell)%>%
  mutate(sample='GSE135851',type='LAM')

LAM1 <- CreateSeuratObject(counts = LAM1.data,
                           meta.data = metadata_t,
                           min.cells = 3,
                           min.features = 200)
LAM1@assays[['RNA']]@layers[['counts']]

features <- LAM1@assays[['RNA']]@features
features <- dimnames(features)[[1]]

NORMAL.data <- Read10X('Normal_Lung/')
hist(log2(NORMAL.data@x),breaks=100)

# 查看log2(x)基因表达情况
par(mfrow=c(1,2))
hist(log2(LAM1.data@x),breaks=100)
hist(log2(NORMAL.data@x),breaks=100)
par(mfrow=c(1,1))

n_gene <- dimnames(NORMAL.data)[[1]]
n_cell <- dimnames(NORMAL.data)[[2]]
metadata_n <- data.frame(row.names = n_cell) %>%
  mutate(sample='GSE135851',type='normallung')
NORMAL <- CreateSeuratObject(counts = NORMAL.data,
                             meta.data = metadata_n,
                             min.cells = 3,
                             min.features = 200)
LAM1_NORMAL <- merge(LAM1,NORMAL)
par(mfrow=c(1,2))
hist(LAM1_NORMAL$nCount_RNA)
hist(LAM1_NORMAL$nFeature_RNA)
par(mfrow=c(1,1))


par(mfrow=c(2,2))
hist(log2(LAM1.data@x),breaks=100, main = 'Figure 1a: log2(gene expression) of LAM')
hist(log2(NORMAL.data@x),breaks=100, main = 'Figure 1b: log2(gene expression) of Normal lung')
hist(LAM1_NORMAL$nCount_RNA, main = 'Figure 1c: Distributions of Gene Counts')
hist(LAM1_NORMAL$nFeature_RNA, main = 'Figure 1d: Distributions of Gene Features')
par(mfrow=c(1,1))
#save(LAM1_NORMAL,file = 'merge_raw_counts.rda')
