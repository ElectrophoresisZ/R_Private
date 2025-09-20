library(Seurat)
library(SingleR)
library(celldex)
library(scrapper)
library(patchwork)
normal <- readRDS('normal_lung_final.rds')
table(normal$seurat_clusters)
normal <- RenameIdents(normal,'0'='type1','1'='type2','2'='type3',
                     '3'='type4','4'='type5','5'='type6',
                     '6'='type7','7'='type8','8'='type9',
                     '9'='type10','10'='type11','11'='type12','12'='type13')
DimPlot(normal, reduction="umap",label=T)
#将Seurat对象normal的默认Assay设置为 'RNA'
DefaultAssay(normal) <- 'RNA'
clusters <- normal$seurat_clusters
humanprimarycell_ref <- celldex::HumanPrimaryCellAtlasData()
#细胞注释
#从Seurat对象normal的默认Assay中提取标准化后的基因表达数据(slot = 'data')
normal_singleR <- GetAssayData(normal,slot = 'data')
pre_annotation <- SingleR(test = normal_singleR,
                          ref = humanprimarycell_ref,
                          labels = humanprimarycell_ref$label.fine,
                          method = 'cluster',
                          clusters = clusters)
table(pre_annotation$labels)
#细胞类型注释的得分热图
plotScoreHeatmap(pre_annotation,clusters = pre_annotation$labels)
new.cluster <- pre_annotation$labels
names(new.cluster) <- levels(normal)
normal <- RenameIdents(normal, new.cluster)
DimPlot(normal,reduction = 'umap',label = T, repel = T)
