Macrophage_cells <- subset(immune.combined, 
                           idents=c('Macrophage:Alveolar','Macrophage:monocyte-derived:IL-4/Dex/TGFb'))
Idents(Macrophage_cells) <- "orig.ident"
table(Idents(Macrophage_cells))
macrophage_gene <- FindMarkers(Macrophage_cells, 
                               ident.1 = "normal_lung",  # 第一组
                               ident.2 = "LAM",  # 第二组
                               assay = "RNA",  
                               logfc.threshold = 0.25, 
                               min.pct = 0.1,  # 基因至少在 10% 的细胞中表达
                               verbose = FALSE)
#save(macrophage_gene,file ='macrophage_gene.Rdata' )
sig_macrophage_gene <- macrophage_gene[macrophage_gene$p_val_adj<0.05 & abs(macrophage_gene$avg_log2FC)>0.5,]
macrophage_gene$sign <- 'no'
macrophage_gene[macrophage_gene$avg_log2FC>0.5,'sign'] <- 'Up'
macrophage_gene[macrophage_gene$avg_log2FC< -0.5,'sign'] <- 'Down'
ggplot(macrophage_gene, aes(x=avg_log2FC,y=-log10(p_val),colour = sign)) +
  geom_point() +
  scale_color_manual(values=c("royalblue4","grey70","firebrick3"))


top_genes <- rownames(sig_macrophage_gene)[1:40]
unique(duplicated(top_genes))

avg.macrophage.cell <- as.data.frame(log1p(AverageExpression(Macrophage_cells, verbose = FALSE)$RNA))
avg.macrophage.cell$gene <- rownames(avg.macrophage.cell)
colnames(avg.macrophage.cell) <- c('LAM','normallung','gene')
p1 <- ggplot(avg.macrophage.cell, aes(normallung,LAM)) + 
  geom_point() + 
  ggtitle("Macrophage cells")
LabelPoints(plot = p1, 
            points = top_genes, 
            repel = TRUE)
FeaturePlot(Macrophage_cells, 
            features = c("CD52","APOE","HBB","CCL2"), 
            split.by = "orig.ident", 
            max.cutoff = 3, #对于归一化数据,表达值超过 3 的会被截断为
            cols = c("grey", "red"))

