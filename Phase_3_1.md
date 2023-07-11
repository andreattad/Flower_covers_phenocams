# Feature selection and models comparison
```r
library(randomForest)
library(crfsuite)
library(varSel)

rm(list=ls())

summary<-data.frame(name_test=c("all","rgb","rgb_vi","rgb_tex","sel"),
                    bands_used=c(rep(NA,5)),
                    f1=c(rep(NA,5))
                    )
# ------------------------------------------------------------
# 0 PREPARE THE TRAINING AND TESTING DATASETS
#------------------------------------------------------------
df<-read.csv(file = "your/folder/path/Phase_2_df_ws_ext_feat/Phase2_df4_ws11.csv",row.names=1 )
df$type<-as.factor(df$type)
df$vari[df$vari==Inf]<-max(df$vari[is.finite(df$vari)])
df$vari[df$vari==-Inf]<-min(df$vari[is.finite(df$vari)])
summary(df$vari)

# SPLIT DATASETS INTO TRAINING AND TESTING 
# SO THAT POINTS OF SAME IMAGE ARE IN THE SAME DATASET
set.seed(1409)
train_index<- sample(seq(1,length(unique(df$imname))),round(length(unique(df$imname))*0.7))
imgst<-unique(df$imname)[train_index]

df$tv<-NULL
train<-df[is.element(df$imname,imgst),c("tv")]<-"t"
test<-df[!is.element(df$imname,imgst),c("tv")]<-"v"
summary(as.factor(df$tv))
train<-df[df$tv=="t",]
test<-df[df$tv=="v",]


#------------------------------------------------------------
# 1. RF MODEL INCLUDING ALL FEATURES
#------------------------------------------------------------
set.seed(1509)  # Setting seed
classifier_RF_ALLB = randomForest(x = train[,c(1:28)],
                             y = as.factor(train$type),
                             ntree = 500,proximity=T)

pred_rf_ALLB = predict(classifier_RF_ALLB, newdata = test)
metrics<-crf_evaluation(pred_rf_ALLB, test$type)
metrics
table(pred_rf_ALLB, test$type)
saveRDS(classifier_RF_ALLB,"your/folder/path/Phase_3_RF_classifiers/RF_all.rds")
 
summary[summary$name_test=="all","bands_used"]<-paste0("'",paste(colnames(train[,c(1:28)]),collapse="','"),"'")
summary[summary$name_test=="all","f1"]<-metrics[["overall"]][["f1_mean"]]
  
  #------------------------------------------------------------
  # 2. RF MODEL INCLUDING RGB BANDS
  #------------------------------------------------------------
 bestBands<-c("R","G","B")
 set.seed(1509)  # Setting seed
 classifier_RF_SELB = randomForest(x = train[,bestBands],
                                   y = as.factor(train$type),
                                   ntree = 500,proximity=T)
 
 pred_rf_SELB = predict(classifier_RF_SELB, newdata = test)
 metrics<-crf_evaluation(pred_rf_SELB, test$type)
 bestBands
 table(pred_rf_SELB, test$type)
 saveRDS(classifier_RF_SELB,"your/folder/path/Phase_3_RF_classifiers/RF_rgb.rds")

 summary[summary$name_test=="rgb","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
 summary[summary$name_test=="rgb","f1"]<-metrics[["overall"]][["f1_mean"]]
 #------------------------------------------------------------
 # 3. RF MODEL INCLUDING RGB BANDS + VEGETATIONAL INDICES
 #------------------------------------------------------------

 bestBands<-c("R","G","B","rgbvi","gli","vari","ngrdi") 
 set.seed(1509)  # Setting seed
 classifier_RF_SELB = randomForest(x = train[,bestBands],
                                   y = as.factor(train$type),
                                   ntree = 500,proximity=T)
 
 pred_rf_SELB = predict(classifier_RF_SELB, newdata = test)
 metrics<- crf_evaluation(pred_rf_SELB, test$type)
 bestBands
table(pred_rf_SELB, test$type)
 saveRDS(classifier_RF_SELB,"your/folder/path/Phase_3_RF_classifiers/RF_rgb_vi.rds")
 
 summary[summary$name_test=="rgb_vi","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
 summary[summary$name_test=="rgb_vi","f1"]<-metrics[["overall"]][["f1_mean"]]
 
 
 #------------------------------------------------------------
 # 4. RF MODEL INCLUDING RGB BANDS + TEXTURE METRICS
 #------------------------------------------------------------
 bestBands<- colnames(df[,c(1:3,8:28)])
 #Build the classifier WITH SELECTED  BANDS
 set.seed(1509)  # Setting seed
 classifier_RF_SELB = randomForest(x = train[,bestBands],
                                   y = as.factor(train$type),
                                   ntree = 500,proximity=T)

 pred_rf_SELB = predict(classifier_RF_SELB, newdata = test)
 metrics<- crf_evaluation(pred_rf_SELB, test$type)
 bestBands
 print(table(pred_rf_SELB, test$type))
 saveRDS(classifier_RF_SELB,"your/folder/path/Phase_3_RF_classifiers/RF_rgb_texm.rds")
 
 summary[summary$name_test=="rgb_tex","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
 summary[summary$name_test=="rgb_tex","f1"]<-metrics[["overall"]][["f1_mean"]]
 
 # ------------------------------------------------------------
 # 5 RF MODEL INCLUDING RGB BANDS + BANDS SELECTED BY SFFS
 #------------------------------------------------------------

 results<-varSelSFFS(
    X = train[,c(1:28)],
    g = as.factor(train$type)
 )
 plot(as.numeric(results$distances),type="b")
 abline(h=sqrt(2),col="red",lty=2)
 abline(v=11,col="green",lty=2)

  bestBands<-names(train)[results$features[11,1:11]]
 #Build the classifier WITH SELECTED  BANDS
 set.seed(1509)  # Setting seed
 classifier_RF_SELB = randomForest(x = train[,bestBands],
                                   y = as.factor(train$type),
                                   ntree = 500,proximity=T)
 
 pred_rf_SELB = predict(classifier_RF_SELB, newdata = test)
 metrics<-crf_evaluation(pred_rf_SELB, test$type)
 bestBands; 
metrics[["bylabel"]][["label"]]
round(metrics[["bylabel"]][["precision"]],2)
metrics[["bylabel"]][["label"]]
round(metrics[["bylabel"]][["recall"]],2)
metrics
 
print(table(pred_rf_SELB, test$type))
 saveRDS(classifier_RF_SELB,"your/folder/path/Phase_3_RF_classifiers/RF_selected.rds")
 plot(as.numeric(results[["distances"]]))
 
 summary[summary$name_test=="sel","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
 summary[summary$name_test=="sel","f1"]<-metrics[["overall"]][["f1_mean"]]

  ###--------------------6---PLOT ALL BANDS
 write.csv(summary,"your/folder/path/Phase_3_RF_classifiers_accuracies.csv")

summary <- summary[order(summary$f1,decreasing = FALSE),]
 barplot(summary$f1,names.arg=summary$name_test,ylim=c(0.70,0.95),xpd=FALSE)
 box(bty="l")
```
