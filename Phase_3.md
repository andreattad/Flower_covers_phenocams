# 3.1 Feature selection and models comparison
We selected a set of best suitable features to optimise processing time, and to reduce redundancy of highly correlated features. 
First, we randomly assigned 70% of images for training, and 30% of images for validation. Then, we compared the accuracies of RF models trained on different subsets of features from the training dataset, including: i) features selected by SFFS, ii) RGB bands alone, iii) RGB bands combined with vegetation indices, iv) RGB bands combined with texture metrics, and v) all features. 

The classifiers were saved in "your/folder/path/Phase_3_RF_classifiers/" and accuracies were saved in "your/folder/path/Phase_3_RF_classifiers_accuracies.csv"

```r
library(randomForest)
library(crfsuite)
library(varSel)
setwd(maindir<-"your/folder/path")

summary<-data.frame(name_test=c("all","rgb","rgb_vi","rgb_tex","sel"),
                    bands_used=c(rep(NA,5)),
                    f1=c(rep(NA,5))
                    )
# ------------------------------------------------------------
# 0 PREPARE THE TRAINING AND TESTING DATASETS
#------------------------------------------------------------
# Load the labelled and sampled points calculated using the best downscaling factor - window size combinations, as resulted from phase 2.2
# Load the dataset of the labelled and sampled pixels calculated using the best df-ws combination (phase 2.2) and replace VARI infinite values
  data<-NA
  data<-read.csv(stringsAsFactors=T,paste0("./Phase_2_df_ws_ext_feat/",files[i]),row.names=1)
  data$vari<-as.numeric(data$vari)
  data$vari[is.na(data$vari)]<-0
  data$vari[data$vari==Inf]<-max(data$vari[is.finite(data$vari)])
  data$vari[data$vari==-Inf]<-min(data$vari[is.finite(data$vari)])

# Split the dataset into training and validation so that points of the same image are in the same subset.
  set.seed(1409)
  train_index<- sample(seq(1,length(unique(data$imname))),round(length(unique(data$imname))*0.7))
  imgst<-unique(data$imname)[train_index]
  data$tv<-NULL #We create a new column that labels each observation as either "training" or "validation" based on the subdataset to which it has been assigned.
  train<-data[is.element(data$imname,imgst),c("tv")]<-"t"
  valid<-data[!is.element(data$imname,imgst),c("tv")]<-"v"
  train<-data[data$tv=="t",]
  valid<-data[data$tv=="v",]
  summary(train$type); summary(valid$type)


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
  saveRDS(classifier_RF_ALLB,"./Phase_3_RF_classifiers/RF_all.rds")
   
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
  saveRDS(classifier_RF_SELB,"./Phase_3_RF_classifiers/RF_rgb.rds")
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
  saveRDS(classifier_RF_SELB,"./Phase_3_RF_classifiers/RF_rgb_vi.rds")
  summary[summary$name_test=="rgb_vi","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
  summary[summary$name_test=="rgb_vi","f1"]<-metrics[["overall"]][["f1_mean"]]
 
#------------------------------------------------------------
# 4. RF MODEL INCLUDING RGB BANDS + TEXTURE METRICS
#------------------------------------------------------------
  bestBands<- colnames(df[,c(1:3,8:28)]) 
  set.seed(1509)  # Setting seed
  classifier_RF_SELB = randomForest(x = train[,bestBands],
                                   y = as.factor(train$type),
                                   ntree = 500,proximity=T)
  pred_rf_SELB = predict(classifier_RF_SELB, newdata = test)
  metrics<- crf_evaluation(pred_rf_SELB, test$type)
  print(table(pred_rf_SELB, test$type))
  saveRDS(classifier_RF_SELB,"./Phase_3_RF_classifiers/RF_rgb_texm.rds")
 
  summary[summary$name_test=="rgb_tex","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
  summary[summary$name_test=="rgb_tex","f1"]<-metrics[["overall"]][["f1_mean"]]
 
# ------------------------------------------------------------
# 5 RF MODEL INCLUDING RGB BANDS + BANDS SELECTED BY SFFS
#------------------------------------------------------------
  # Select the most informative bands using varSelSFFS
    results<-varSelSFFS(X = train[,c(1:28)],g = as.factor(train$type))
    plot(as.numeric(results$distances),type="b")
    abline(h=sqrt(2),col="red",lty=2)
    abline(v=11,col="green",lty=2)
    bestBands<-names(train)[results$features[11,1:11]]

  # Build the classifier WITH SELECTED  BANDS
    set.seed(1509)  # Setting seed
    classifier_RF_SELB = randomForest(x = train[,bestBands],
                                   y = as.factor(train$type),
                                   ntree = 500,proximity=T)
 
    pred_rf_SELB = predict(classifier_RF_SELB, newdata = test)

  # Print the results
    metrics<-crf_evaluation(pred_rf_SELB, test$type)
    metrics                  
    print(table(pred_rf_SELB, test$type))
    saveRDS(classifier_RF_SELB,"./Phase_3_RF_classifiers/RF_selected.rds")
    plot(as.numeric(results[["distances"]]))
    summary[summary$name_test=="sel","bands_used"]<-paste0("'",paste( bestBands,collapse="','"),"'")
    summary[summary$name_test=="sel","f1"]<-metrics[["overall"]][["f1_mean"]]

###--------------------Save the accuracies summary
 write.csv(summary,"./Phase_3_RF_classifiers_accuracies.csv")


```

