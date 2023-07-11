# Model comparison
Here, we tested the influence of window size and downscaling factor on classification accuracy.
We randomly assigned pixels from 70% of the images for training, and 30% for validation for each DF-WS combination, then we trained a random forest (RF) classifier using the “randomForest” function of the “randomForest” package (Liaw & Wiener, 2002). Accuracies were reported and saves as .csv file. 
The metric to calculate the accuracy of the RF classifiers was the mean F1 score of the five classes. The F1 score is derived from precision and recall metrics. The precision is intuitively the ability of the classifier not to label a sampled pixel as positive when it is negative, whereas the recall is the ability of the classifier to find all the positive sampled pixels. Metrics equations are reported in the main text. 
```
library(randomForest)
library(crfsuite)
library(bigstatsr)
library(readr)
library(dplyr)
library(raster)
rm(list=ls())

files<-list.files(path="your/folder/path/Phase_2_df_ws_ext_feat/",
                  pattern="Phase2")
res<-data.frame(f1_mean=c(rep(NA,length(files))),
                df=c(rep(NA,length(files))),
                ws=c(rep(NA,length(files))),
                npoints=c(rep(NA,length(files))))

for (i in 1:length(files) ){
   print(i);   print(files[i])
  # prepare train e test samples
  data<-NA
  data<-read.csv(stringsAsFactors=T,paste0("your/folder/path/Phase_2_df_ws_ext_feat/",files[i]),row.names=1)
  data$vari<-as.numeric(data$vari)
  data$vari[is.na(data$vari)]<-0
  data$vari[data$vari==Inf]<-max(data$vari[is.finite(data$vari)])
  data$vari[data$vari==-Inf]<-min(data$vari[is.finite(data$vari)])
  #split the dataset into training and testing so that points of the same imageare in the same subset.
    set.seed(1409)
    train_index<- sample(seq(1,length(unique(data$imname))),round(length(unique(data$imname))*0.7))
    imgst<-unique(data$imname)[train_index]
    data$tv<-NULL
    train<-data[is.element(data$imname,imgst),c("tv")]<-"t"
    test<-data[!is.element(data$imname,imgst),c("tv")]<-"v"
    train<-data[data$tv=="t",]
    test<-data[data$tv=="v",]
    summary(train$type); summary(test$type)

  #Build the classifier
    set.seed(1509)  # Setting seed
  classifier_RF = randomForest(x = train[,c(1:28)],
                               y = as.factor(train$type),
                               ntree = 500,proximity=T)
  pred_rf = predict(classifier_RF, newdata = test)
  # pred_rf[pred_rf=="Soil"]<-"Green_vegetation"
  # test$type[test$type=="Soil"]<-"Green_vegetation"
  # test$type<-factor(test$type)
  print(crf_evaluation(pred_rf, test$type))
  print(table(pred_rf, test$type))
  res$f1_mean[i]<-crf_evaluation(pred_rf, test$type)[["overall"]][["f1_mean"]]
  res$df[i]<-parse_number(substr(files[i],10,11))
  res$ws[i]<-parse_number(substr(files[i],14,15))
  res$npoints[i]<-nrow(data)

}
write.csv(res,file="your/folder/path/Phase_2_dfws_accuracy.csv")
```
