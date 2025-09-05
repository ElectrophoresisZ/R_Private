# 数据准备
rm(list = ls())
load('GSE163151_updated.Rdata')
load('perm.pVal.Rdata')
load('Original_Trainset.Rdata')
load('Testset.Rdata')
load('Train_subset.Rdata')
load('Validation_set.Rdata')

#perm.pVal.updated <- perm.pVal[perm.pVal$`COVID-19.pval`<0.25,]
selected_genes <- row.names(perm.pVal)
train_GSE163151 <- train_GSE163151[,c('Response', selected_genes)]
test_GSE163151 <- test_GSE163151[,c('Response', selected_genes)]
train_subset_GSE163151 <- train_subset_GSE163151[,c('Response', selected_genes)]
validation_GSE163151 <- validation_GSE163151[,c('Response', selected_genes)]

# 随机森林建模
library(randomForest)

train.x <- train_subset_GSE163151[,-which(names(train_subset_GSE163151) == 'Response')]
train.y <- train_subset_GSE163151$Response
rf.model <- randomForest(train.x, train.y, importance = T, ntree = 500)
importance(rf.model)
rf.model
#save(rf.model, file = 'rf_model_optimized.Rdata')

valid.x <- validation_GSE163151[,-which(names(validation_GSE163151) == 'Response')]
valid.y <- validation_GSE163151$Response

# ROC曲线
library(ROCR)
pred.rf.prob <- predict(rf.model, valid.x, type = 'prob')
pred.obj.rf <- prediction(pred.rf.prob[,2], valid.y)
perf.auc.rf <- performance(pred.obj.rf, 'auc')
auc.rf <- perf.auc.rf@y.values[[1]]
perf.roc.rf <- performance(pred.obj.rf, 'tpr','fpr')
plot(perf.roc.rf, col = 'blue', lwd = 2)
abline(0,1,lty=2)
legend('bottomright', legend = paste('AUC: ', round(auc.rf, 4)), 
       col = 'tomato', lwd=2, bty = 'n') # 不绘制任何边框

#约登指数
tpr <- perf.roc.rf@y.values[[1]]
fpr <- perf.roc.rf@x.values[[1]]
threshold <- perf.roc.rf@alpha.values
j <- tpr - fpr
threshold[[1]][which.max(j)]

# 混淆矩阵
pred.rf.repsonse <- predict(rf.model, newdata = valid.x)
confumatrix.2 <-caret::confusionMatrix(valid.y, pred.rf.repsonse)
confumatrix.2
#save(confumatrix.2, file = 'confusionmatrix2.Rdata')


# 交叉验证
library(caret)
ctrl <- trainControl(method = 'cv', number = 5,)

library(rpart)
rf.cv <- train(Response ~.,
               data = train_GSE163151,
               method = 'rf',
               trControl = ctrl)
rf.cv
rf.cv$resample
rf.cv$finalModel

#save(rf.cv, file = 'rf_cv_model.Rdata')
load('rf_cv_model.Rdata')

# Kappa>0.6即为良好，>0.8为极好
# P-Value [Acc > NIR]:检验模型的准确性是否显著高于无信息率, P值远小于0.05，说明模型的性能提升在统计上是极其显著的
pred.cv <- predict(rf.cv, newdata = test_GSE163151[,-1])
confumatrix.3 <- caret::confusionMatrix(test_GSE163151$Response, pred.cv)
confumatrix.3
#save(confumatrix.3, file = 'confusionmatrix3.Rdata')

pred.cv.train.prob <- predict(rf.cv, newdata = train_GSE163151[,-1], type = 'prob')
pred.cv.prob <- predict(rf.cv, newdata = test_GSE163151[,-1], type = 'prob')
library(ggplot2)
theme_set(theme_bw() +
            theme(axis.text.x = element_text(size = 10,  color = "black"),
                  axis.text.y = element_text(size = 10, color = "black"),
                  axis.title.x = element_text(size = 12, face = "bold",color = "black"),
                  axis.title.y = element_text(size = 12, face = "bold", color = "black"),
                  #panel.grid.major = element_blank(),
                  #panel.grid.minor = element_blank(),
                  panel.border = element_rect(color = "black", linewidth = 1),  # Set border size here
                  #legend.position = "none",
                  plot.title = element_text(hjust = 0.5, size = 12),  
                  strip.text = element_text(size = 12, face = "bold"),  
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 10)))
library(pROC)
rocc.train <- roc(train_GSE163151$Response, pred.cv.train.prob[,1])
plot(rocc.train, 
     print.auc=TRUE, 
     auc.polygon=TRUE,  # 用颜色填充 ROC 曲线下方的区域
     mar = c(5,4,5,4), # 调整边距 [下,左,上,右]
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue", 
     grid=c(0.1, 0.1),   # 在图中添加网格线，以辅助读图。
     grid.col=c("green", "red"), 
     print.thres=T, 
     main ='ROC of Training Set',
     cex.main = 1.5,
     cex.lab = 1.2) 
rocc <- roc(test_GSE163151$Response, pred.cv.prob[,1])
plot(
  rocc,
  col = "#2E86C1",              # ROC 曲线颜色（蓝色调）
  lwd = 2.5,                    # 线条宽度
  main = "ROC Curve of Testing Set", # 标题
  cex.main = 2,               # 标题字体大小
  cex.lab = 1.5,                # 坐标轴标题字体大小
  cex.axis = 1.0,               # 坐标轴刻度字体大小
  print.thres=T, 
  legacy.axes = TRUE,           # 坐标轴方向改成 FPR vs TPR（更常见）
  print.auc = TRUE,             # 打印 AUC
  print.auc.cex = 1.2,          # AUC 字体大小
  print.auc.col = "black",      # AUC 字体颜色
  grid = FALSE                  # 不加花哨网格
)

# 添加对角线 (随机分类器参考线)
#abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 1.5)

# 添加 AUC 阴影区域（更淡雅）
plot(
  rocc,
  col = "#2E86C1",
  lwd = 2.5,
  add = TRUE,
  print.auc = FALSE,
  auc.polygon = TRUE,
  auc.polygon.col = rgb(46,134,193,50, maxColorValue=255) # 半透明蓝色
)

# 校准曲线
library(ggplot2)
library(dplyr)

# 在每个点附近做局部加权回归，把这些局部拟合结果拼接起来，得到一条平滑曲线
loess.response <- ifelse(GSE163151_updated$Response == 'COVID-19', 1, 0)
all.pred.cv.prob <- predict(rf.cv, newdata = GSE163151_updated[,-1], type = 'prob')
loess.data <- data.frame(pred = all.pred.cv.prob[,1], obs = loess.response)
# p-value > 0.05：无法拒绝原假设，没有足够证据表明模型校准不佳。这是一个“好”结果
library(ResourceSelection)
h.test <- hoslem.test(loess.response, all.pred.cv.prob[,1], g=5)
p_value <- h.test$p.value

ideal_line <- data.frame(x = c(0, 1), y = c(0, 1))
ggplot(loess.data, aes(pred, obs)) +
  geom_point(alpha = 0.5, color = 'grey50') +
  geom_smooth(aes(color = "loess fit"), method = 'loess', se = TRUE) +
  geom_line(data = ideal_line, aes(x = x, y = y, color = "ideal"), linetype = "dashed") +
  xlim(0, 1) + ylim(0, 1) +
  scale_color_manual(
    name = 'Legend',
    values = c("loess fit" = "#2255AD", "ideal" = "#B82613"),
    labels = c("LOESS fit", "Ideal")
  ) +
  labs(x = "Predicted Probability", y = "Observed Frequency",
       title = "Calibration Plot (LOESS Smoothed)") +
  theme_bw() +
  theme(legend.position =  c(0.8, 0.2),
        legend.text = element_text(size = 10, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(hjust = 0.5, size = 16, face = "bold")
        ) +
  annotate("text",
           x = 0.1, y = 0.9, # 将文本放在图的左上角
           label = paste("H-L Test p-value =", round(p_value,3)),
           size = 7, color = "#268C40")


#--------------------------------------------------------------------------------------------

# 置换检验
load('rf.perm.Rdata')
summary(rf.perm)
importance(rf.perm)
plotImportance(rf.perm)

imp <- importance(rf.perm, type = 1) # type=1 对于MeanDecreaseAccuracy
imp_df <- data.frame(
  Feature = rownames(imp),
  Importance = imp[,5],
  Pvalue = imp[,6]
)
imp_df <- imp_df %>%
  mutate(Significant = ifelse(Pvalue < 0.05, "p<0.05", "p>=0.05"))

library(ggplot2)
ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance, fill = Significant)) + #按重要性排序
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("p<0.05" = "tomato", "p>=0.05" = "skyblue")) +
  theme_bw() +
  labs(
    title = "Permutation Importance by P-value",
    x = "Feature",
    y = "Importance",
    fill = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    legend.text = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.y = element_text(hjust = 0.5, size = 16, face = "bold")
  )







