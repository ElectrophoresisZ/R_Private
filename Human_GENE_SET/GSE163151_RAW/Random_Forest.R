# 数据准备
rm(list=ls())
selected_gene <- read.table('GSE163151_select.txt')
load('40selected_genes.Rdata')
load('pdata_GSE163151_subset.Rdata')
GSE163151_plot_trans <- t(GSE163151_plot)
GSE163151_plot_trans <- scale(GSE163151_plot_trans)
Y <- data.frame(Response = pdata_GSE163151_subset$`disease state:ch1`)
row.names(Y) <- rownames(GSE163151_plot_trans)
GSE163151_plot_trans <- cbind(Y, GSE163151_plot_trans)
GSE163151_plot_trans$Response <- as.factor(GSE163151_plot_trans$Response)

# 分层抽样
set.seed(123)
library(caret)
trainindex.ori <- createDataPartition(GSE163151_plot_trans$Response, p=0.8, list = F)
train_GSE163151 <- GSE163151_plot_trans[trainindex.ori,]
test_GSE163151 <- GSE163151_plot_trans[-trainindex.ori,]

trainindex <- createDataPartition(train_GSE163151$Response, p=0.8, list = F)
train_subset_GSE163151 <- train_GSE163151[trainindex,]
validation_GSE163151 <- train_GSE163151[-trainindex,]

save(train_GSE163151, file = 'Original_Trainset.Rdata')
save(test_GSE163151, file = 'Testset.Rdata')
save(train_subset_GSE163151, file = 'Train_subset.Rdata')
save(validation_GSE163151, file = 'Validation_set.Rdata')

# 随机森林建模
library(randomForest)

train.x <- train_subset_GSE163151[,-which(names(train_subset_GSE163151) == 'Response')]
train.y <- train_subset_GSE163151$Response
rf.model <- randomForest(train.x, train.y, importance = T, ntree = 500)
importance(rf.model)
rf.model
#save(rf.model, file = 'rf.model.Rdata')

valid.x <- validation_GSE163151[,-which(names(validation_GSE163151) == 'Response')]
valid.y <- validation_GSE163151$Response

# ROC曲线
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
confumatrix.1 <-caret::confusionMatrix(valid.y, pred.rf.repsonse)
#save(confumatrix.1, file = 'confusionmatrix1.Rdata')

# 置换检验
library(rfPermute)
set.seed(123)
rf.perm <- rfPermute(Response~.,
                     data = train_subset_GSE163151,
                     ntree = 500,
                     nperm = 100)
#save(rf.perm, file = 'rf.perm.Rdata')
rf.perm$rf
summary(rf.perm)
perm.matrix <- importance(rf.perm)
perm.dataframe <- as.data.frame(perm.matrix)
perm.pVal <- perm.dataframe[which(perm.dataframe$MeanDecreaseAccuracy.pval < 0.05),]
#save(perm.pVal, file = 'perm.pVal.Rdata')

selected_genes <- row.names(perm.pVal)
GSE163151_updated <- GSE163151_plot_trans[,selected_genes]
GSE163151_updated <- cbind(Y, GSE163151_updated)
GSE163151_updated$Response <- as.factor(GSE163151_updated$Response)
#save(GSE163151_updated, file = 'GSE163151_updated.Rdata')
