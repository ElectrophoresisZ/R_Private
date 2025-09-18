pbmc <- readRDS('pbmc_LAM1_final.rds')
library(SingleR)
library(celldex)
library(scrapper)
library(cowplot)
library(ggplot2)
library(readxl)
cellmarker <- read_xlsx('Cell_marker_Human.xlsx')

res <- celldex::HumanPrimaryCellAtlasData()
pbmc_for_SingleR <- GetAssayData(pbmc, slot="data")
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = res, labels = res$label.fine)
pbmc.hesc
visual_annoatation <- as.data.frame(table(pbmc.hesc$labels,pbmc@meta.data$seurat_clusters))
visual_annoatation[which(visual_annoatation$Freq>30),]

#-------------------------------------------------------------------------------
table(pbmc$seurat_clusters)
pbmc <- RenameIdents(pbmc,'0'='type1','1'='type2','2'='type3',
                          '3'='type4','4'='type5','5'='type6',
                          '6'='type7','7'='type8','8'='type9',
                          '9'='type10','10'='type11','11'='type12','12'='type13')
DimPlot(pbmc,reduction = "umap",label = T)
DefaultAssay(pbmc) <- 'RNA'
pbmc

pbmc_singler <- GetAssayData(pbmc,slot = "data") # 提取表达矩阵
clusters <- pbmc$seurat_clusters # 获取现在的细胞类群名字
table(clusters)
pred.immune <- SingleR(test = pbmc_singler, #你的表达矩阵
                       ref = res, # 你的注释文件
                       labels = res$label.fine,
                       #因为样本主要为免疫细胞（而不是全部细胞），因此设置为label.fine
                       method = "cluster", 
                       clusters = clusters)
table(pred.immune$labels)
plotScoreHeatmap(pred.immune, clusters = pred.immune$labels)
new.clusterID <- pred.immune$labels
names(new.clusterID) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.clusterID)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE)
