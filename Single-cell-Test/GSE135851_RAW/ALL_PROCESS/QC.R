library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)


load('merge_raw_counts.rda')
#直方图
LAM1_NORMAL$Percent.Mito <- PercentageFeatureSet(LAM1_NORMAL,pattern = '^MT-')
hist(LAM1_NORMAL$Percent.Mito)

LAM1_NORMAL$Percent.ERCC <- PercentageFeatureSet(LAM1_NORMAL,pattern = '^ERCC')
hist(LAM1_NORMAL$Percent.ERCC)

LAM1_NORMAL$Percent.Ribo <- PercentageFeatureSet(LAM1_NORMAL,pattern = '^RP[SL]')
hist(LAM1_NORMAL$Percent.Ribo)


par(mfrow = c(1,3))
# Mito&ERCC&Ribo
hist(LAM1_NORMAL$Percent.Mito, 
     main = "Figure2a: Mitochondrial Gene Percentage",
     xlab = "Percentage of Genes", 
     ylab = "Number of Cells")
hist(LAM1_NORMAL$Percent.ERCC,main = 'Figure2b: External Gene Percentage',
     xlab = "Percentage of Genes", 
     ylab = "Number of Cells")
hist(LAM1_NORMAL$Percent.Ribo,main = 'Figure2c: Ribosomal Gene Percentage',
     xlab = "Percentage of Genes", 
     ylab = "Number of Cells")


library(cowplot)
plot_grid(p1, p2, p3, p4,
          nrow = 2, ncol = 2,
          labels = c("2d", "2e", "2f", "2g"))


#小提琴图
p1 <- VlnPlot(LAM1_NORMAL,features = 'nFeature_RNA',group.by = 'type',
        alpha = 0.1,pt.size = 0.01) + # 控制每个数据点的大小
  geom_hline(yintercept = 2000,color = 'red')+  # 添加一条水平线;指定水平线的y轴位置为2000
  scale_fill_igv() # 提供了一套预设的颜色方案

p2 <- VlnPlot(LAM1_NORMAL,features = 'nCount_RNA',group.by = 'type',
        alpha = 0.1,pt.size = 0.01)+
  geom_hline(yintercept = 5000,color = 'red')+
  scale_fill_igv()

p3 <- VlnPlot(LAM1_NORMAL,features = 'Percent.Ribo',group.by = 'type',
        alpha = 0.1,pt.size = 0.01)+
  geom_hline(yintercept = 2000,color = 'red')+
  scale_fill_igv()

p4 <- VlnPlot(LAM1_NORMAL,features = 'Percent.Mito',group.by = 'type',
        alpha = 0.1,pt.size = 0.01)+
  geom_hline(yintercept = 2000,color = 'red')+
  scale_fill_igv()

# 细胞点图
# UMI 是 Unique Molecular Identifier（唯一分子标识符）的缩写。
correlation <- cor(LAM1_NORMAL$nFeature_RNA, LAM1_NORMAL$nCount_RNA)
cor_text <- paste0("R = ", round(correlation, 3))

FeatureScatter(object = LAM1_NORMAL,
                     group.by = 'type',
                     raster = F, #是否将点栅格化
                     shuffle = T, #是否随机打乱点的绘制顺序
                     pt.size = 0.05,
                     feature1 = "nFeature_RNA",
                     feature2 = "nCount_RNA")+
  scale_color_igv()+
  guides(color = guide_legend(override.aes = list(size = 4))) + #让图例中的点变大
  labs(title = "2h: Feature Scatter Plot (Genes & UMI counts)",
       subtitle = paste0("Correlation Coefficient: ", cor_text)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"))

correlation <- cor(LAM1_NORMAL$Percent.Mito, LAM1_NORMAL$Percent.Ribo)
cor_text <- paste0("R = ", round(correlation, 3))

FeatureScatter(object = LAM1_NORMAL,
                     group.by = 'type',
                     raster = F,
                     shuffle = T,
                     pt.size = 0.05,
                     feature1 = "Percent.Ribo",
                     feature2 = "Percent.Mito")+
  scale_color_igv()+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "2i: Feature Scatter Plot (Ribo & UMI Mito)",
       subtitle = paste0("Correlation Coefficient: ", cor_text)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"))
  

LAM1_NORMAL_scaled <- subset(LAM1_NORMAL,nCount_RNA>1000)
LAM1_NORMAL_scaled <- subset(LAM1_NORMAL_scaled,nFeature_RNA>500)
LAM1_NORMAL_scaled = subset(LAM1_NORMAL_scaled,Percent.Mito<20)
LAM1_NORMAL_scaled = subset(LAM1_NORMAL_scaled,Percent.Ribo<40)
dim(LAM1_NORMAL_scaled)
#save(LAM1_NORMAL_scaled,file = 'LAM_NORMAL_filtered.rda')
