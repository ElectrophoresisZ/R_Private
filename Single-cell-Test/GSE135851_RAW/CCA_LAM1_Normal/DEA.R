rm(list = ls())
library(Seurat)
library(SeuratData)
library(patchwork)
library(metap)

LAM1 <- readRDS('LAM1_final.rds')
LAM1@meta.data$orig.ident <- 'LAM'
normal <- readRDS('normal_lung_final.rds')
#检查基因名是否一致
identical(rownames(LAM1), rownames(normal))
#检查细胞数量
ncol(LAM1);ncol(normal)
LAM.list <- list(LAM=LAM1,Health=normal)
LAM.list <- lapply(X=LAM.list,FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method='vst',nfeatures=2000)
  })
#选择相同的高变基因
features <- SelectIntegrationFeatures(object.list = LAM.list)
#找锚点
#锚点是细胞对（来自不同数据集的细胞），它们在生物学上相似（例如，表达模式相近），可以用来对齐数据集
immune.anchors <- FindIntegrationAnchors(object.list = LAM.list, 
                                         anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
immune.combined

#下游分析
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#指定计算的主成分（PCs）数量为 30
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)#关闭详细输出
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
#saveRDS(immune.combined,file = 'LAM1_NORMAL_LUNG.rds')

# 降维结果
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident") #根据元数据中的 "stim" 列进行着色或分组
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE) #避免标签重叠
p1 + p2

DefaultAssay(immune.combined) <- "RNA"
# 合并数据层
immune.combined <- JoinLayers(immune.combined, assay = "RNA")
nk.markers <- FindConservedMarkers(immune.combined, 
                                   ident.1 = "Epithelial_cells:bronchial", 
                                   grouping.var = "orig.ident", 
                                   verbose = FALSE)
Epithelial_marker <- nk.markers
#save(Epithelial_marker,file = 'Epithelial_marker.Rdata')
head(Epithelial_marker)
#前10%的最高表达值以下全部视为低表达并设为统一颜色
FeaturePlot(immune.combined, 
            features = rownames(head(Epithelial_marker,9)), 
            min.cutoff = "q9") 
DotPlot(immune.combined, features = rownames(head(Epithelial_marker,9)), 
        cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()


#-----------------------------------------------------------------------------
#注释
library(SingleR)
library(celldex)
library(scrapper)
table(immune.combined$seurat_clusters)
humanprimarycell.ref <- HumanPrimaryCellAtlasData()
clusters <- immune.combined$seurat_clusters
integrated_singler <- GetAssayData(immune.combined, slot = 'data')
integrated_hesc <- SingleR(test = integrated_singler,
                           ref = humanprimarycell.ref,
                           labels = humanprimarycell.ref$label.fine,
                           method = 'cluster',
                           clusters = clusters)
#save(integrated_hesc, file = 'singleR_result.Rdata')
table(integrated_hesc$labels)
plotScoreHeatmap(integrated_hesc, clusters = integrated_hesc$labels)
new.cluster <- integrated_hesc$labels
names(new.cluster) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined,new.cluster)
DimPlot(immune.combined, reduction = 'umap',label = T, repel = T)

immune.combined@meta.data$singleR_labels <- Idents(immune.combined)
#-----------------------------------------------------------------------------
#DEA
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

epithelial.cells <- subset(immune.combined, idents='Epithelial_cells:bronchial')
Idents(epithelial.cells) <- "orig.ident"
#计算每个基因在不同细胞群体中的平均 RNA 表达量
avg.epithelial.cells <- as.data.frame(log1p(AverageExpression(epithelial.cells, verbose = FALSE)$RNA)) 
avg.epithelial.cells$gene <- rownames(avg.epithelial.cells)
colnames(avg.epithelial.cells) <- c('LAM','normallung','gene')
p1 <- ggplot(avg.epithelial.cells, aes(normallung,LAM)) + 
  geom_point() + 
  ggtitle("Epithelial cells")
LabelPoints(plot = p1, 
            points = top_genes, 
            repel = TRUE)

FeaturePlot(epithelial.cells, 
            features = c("SFTPC", "SFTPD", "HBB","MT-ND2","NME1-NME2"), 
            split.by = "orig.ident", 
            max.cutoff = 3, #对于归一化数据,表达值超过 3 的会被截断为
            cols = c("grey", "red"))








#-------------------------------------------------------------------------------
table(Idents(epithelial.cells))
# orig.ident有两个组：Control 和 Treatment
de_genes <- FindMarkers(epithelial.cells, 
                        ident.1 = "normal_lung",  # 第一组
                        ident.2 = "LAM",  # 第二组
                        assay = "RNA",  
                        logfc.threshold = 0.25, 
                        min.pct = 0.1,  # 基因至少在 10% 的细胞中表达
                        verbose = FALSE)
epithelial_gene <- de_genes
#save(epithelial_gene,file = 'epithelial_gene.Rdata')
head(de_genes)
sig_de_genes <- de_genes[de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) > 0.5, ]

ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot")

library(pheatmap)
# 筛选显著差异基因
top_genes <- rownames(sig_de_genes)[1:40] 
unique(duplicated(top_genes))
expr_matrix <- GetAssayData(epithelial.cells, assay = "RNA", slot = "data")
expr_matrix <- as.data.frame(expr_matrix)
expr_matrix_select <- expr_matrix[top_genes,]
annotation_col <- data.frame(Type=epithelial.cells@meta.data$orig.ident)
rownames(annotation_col) <- colnames(expr_matrix_select)
pheatmap(expr_matrix_select, scale = 'row',
         show_rownames = T, show_colnames = F,
         cluster_rows = F, cluster_cols = T,
         annotation = annotation_col,
         color = colorRampPalette(c("royalblue4","grey","firebrick3"))(500))


