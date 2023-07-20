# 3.1 Feature selection and model comparison
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
# PREPARE THE TRAINING AND TESTING DATASETS
#------------------------------------------------------------
# Load the labelled and sampled points calculated using the best df-ws combinations (phase 2.2)
  data<-NA
  data<-read.csv(stringsAsFactors=T,paste0("./Phase_2_df_ws_ext_feat/",files[i]),row.names=1)
  data$vari<-as.numeric(data$vari)
  data$vari[is.na(data$vari)]<-0
  data$vari[data$vari==Inf]<-max(data$vari[is.finite(data$vari)])
  data$vari[data$vari==-Inf]<-min(data$vari[is.finite(data$vari)])

# Split the dataset into training and validation so that points of the same image are in the same subset
  set.seed(1409)
  train_index<- sample(seq(1,length(unique(data$imname))),round(length(unique(data$imname))*0.7))
  imgst<-unique(data$imname)[train_index]
  data$tv<-NULL
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
   
  summary[summary$name_test=="all","bands_used"]<-paste0("'",
                                        paste(colnames(train[,c(1:28)]),collapse="','"),"'")
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

# 3.2 Processing time calculation
We calculated the processing time for the calculation of all features and subsequent image classification for one image.
The durations were reported in "your/folder/path/Phase_3_RF_classifiers_accuracies.csv". The feature combination providing the best trade-off between accuracy and processing time was selected for the RF final classifier compilation. 

```r
library(randomForest)
library(glcm)
library(raster)
library(terra)
library(crfsuite)
setwd(maindir<-"your/folder/path")

summary<-read.csv("./Phase_3_RF_classifiers_accuracies.csv",row.names=1)
summary$durations<-NA

df<-4
ws<-11
load(paste("./ROISREFS_2014/A/010/ROI/roi.data.Rdata",sep=""))


#------------------------------------------------------------
#  1 PROCESSING TIME OF RF MODELS INCLUDING ALL BANDS
#------------------------------------------------------------
  classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
  classifier_RF<-readRDS("./Phase_3_RF_classifiers/RF_all.rds")
  time0 <- proc.time()[3]
    imgr<-brick("./IMGS_2014/010/SiteJEA010_201404121334.jpg")
    imgr <- terra::aggregate(imgr, df)
    imgr[[4]]<-((imgr[[2]]*imgr[[2]])-(imgr[[1]]*imgr[[3]]))/((imgr[[2]]*imgr[[2]])+(imgr[[1]]*imgr[[3]]))#RGBVI
    imgr[[5]]<-((2*imgr[[2]])-imgr[[1]]-imgr[[3]])/((2*imgr[[2]])+imgr[[1]]+imgr[[3]])#GLI
    imgr[[6]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]]-imgr[[3]])#VARI
    imgr[[7]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]])#NGRDI
    imgr[[6]][imgr[[6]]==Inf]<-max(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))]) #to avoid +inf value in VARI
    imgr[[6]][imgr[[6]]==-Inf]<-min(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))])
    glcm.red <- glcm(imgr[[1]],
                     window = c(ws, ws),#deve essere dispari
                     na_opt="center",min_x=0,max_x=255,
                     shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                     statistics = c("variance","homogeneity","contrast","entropy",
                                    "dissimilarity", "second_moment","mean"))
    glcm.green <- glcm(imgr[[2]],
                     window = c(ws, ws),#deve essere dispari
                     na_opt="center",min_x=0,max_x=255,
                     shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                     statistics = c("variance","homogeneity","contrast","entropy",
                                    "dissimilarity", "second_moment","mean"))
    glcm.blue <- glcm(imgr[[3]],
                     window = c(ws, ws),#deve essere dispari
                     na_opt="center",min_x=0,max_x=255,
                     shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                     statistics = c("variance","homogeneity","contrast","entropy",
                                    "dissimilarity", "second_moment","mean"))
    imgm<-mask(stack(imgr,glcm.red,glcm.green,glcm.blue),roi.data[[1]]$polygons)
    names(imgm)<-c("R","G","B","rgbvi","gli","vari","ngrdi",
           "Rvariance","Rhomogeneity","Rcontrast","Rentropy","Rdissimilarity","Rsecond_moment","Rmean",
           "Gvariance","Ghomogeneity","Gcontrast","Gentropy","Gdissimilarity","Gsecond_moment","Gmean",
           "Bvariance","Bhomogeneity","Bcontrast","Bentropy","Bdissimilarity","Bsecond_moment","Bmean")
    classified<-terra::predict(imgm,classifier_RF,na.rm=T)
    time1 <- proc.time()[3]
    summary$durations[1]<-duration <- as.numeric(time1-time0)
    print(duration)
      
#------------------------------------------------------------
#  2 PROCESSING TIME OF RF MODELS INCLUDING RGB BANDS
#------------------------------------------------------------
  classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
  classifier_RF<-readRDS("./Phase_3_RF_classifiers/RF_rgb.rds")
  time0 <- proc.time()[3]
    imgr<-brick("./IMGS_2014/010/SiteJEA010_201404121334.jpg")
    imgr <- terra::aggregate(imgr, df)
    imgm<-mask(imgr,roi.data[[1]]$polygons)
    names(imgm)<-c("R","G","B")
    classified<-terra::predict(imgm,classifier_RF,na.rm=T)
    time1 <- proc.time()[3]
    summary$durations[2]<-duration <-  as.numeric(time1-time0)
    print(duration) 

#------------------------------------------------------------
#  3 PROCESSING TIME OF RF MODELS INCLUDING RGB+VI
#------------------------------------------------------------
  classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
  classifier_RF<-readRDS("./Phase_3_RF_classifiers/RF_rgb_vi.rds")
  time0 <- proc.time()[3]
    imgr<-brick("./IMGS_2014/010/SiteJEA010_201404121334.jpg")
    imgr <- terra::aggregate(imgr, df)
    imgr[[4]]<-((imgr[[2]]*imgr[[2]])-(imgr[[1]]*imgr[[3]]))/((imgr[[2]]*imgr[[2]])+(imgr[[1]]*imgr[[3]]))#RGBVI
    imgr[[5]]<-((2*imgr[[2]])-imgr[[1]]-imgr[[3]])/((2*imgr[[2]])+imgr[[1]]+imgr[[3]])#GLI
    imgr[[6]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]]-imgr[[3]])#VARI
    imgr[[7]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]])#NGRDI
    imgr[[6]][imgr[[6]]==Inf]<-max(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))]) #to avoid +inf value in VARI
    imgr[[6]][imgr[[6]]==-Inf]<-min(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))])
    imgm<-mask(imgr,roi.data[[1]]$polygons)
    names(imgm)<-c("R","G","B","rgbvi","gli","vari","ngrdi")
    classified<-terra::predict(imgm,classifier_RF,na.rm=T)
    time1 <- proc.time()[3]
    summary$durations[3]<-duration <- as.numeric(time1-time0)
    print(duration)

#------------------------------------------------------------
#  4 PROCESSING TIME OF RF MODELS INCLUDING RGB+TEXTURE METRICS
#------------------------------------------------------------
  classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
  classifier_RF<-readRDS("./Phase_3_RF_classifiers/RF_rgb_texm.rds")
  time0 <- proc.time()[3]
    imgr<-brick("./IMGS_2014/010/SiteJEA010_201404121334.jpg")
    imgr <- terra::aggregate(imgr, df)
    glcm.red <- glcm(imgr[[1]],
                     window = c(ws, ws),#deve essere dispari
                     na_opt="center",min_x=0,max_x=255,
                     shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                     statistics = c("variance",       "homogeneity","contrast",        "entropy",
                                    "dissimilarity", "second_moment","mean"))
    glcm.green <- glcm(imgr[[2]],
                       window = c(ws, ws),#deve essere dispari
                       na_opt="center",min_x=0,max_x=255,
                       shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                       statistics = c("variance",       "homogeneity","contrast",        "entropy",
                                      "dissimilarity", "second_moment","mean"))
    glcm.blue <- glcm(imgr[[3]],
                      window = c(ws, ws),#deve essere dispari
                      na_opt="center",min_x=0,max_x=255,
                      shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                      statistics = c("variance",       "homogeneity","contrast",        "entropy",
                                     "dissimilarity", "second_moment","mean"))
    imgm<-mask(stack(imgr,glcm.red,glcm.green,glcm.blue),roi.data[[1]]$polygons)
    names(imgm)<-c("R","G","B",
           "Rvariance","Rhomogeneity","Rcontrast","Rentropy","Rdissimilarity","Rsecond_moment","Rmean",
           "Gvariance","Ghomogeneity","Gcontrast","Gentropy","Gdissimilarity","Gsecond_moment","Gmean",
           "Bvariance","Bhomogeneity","Bcontrast","Bentropy","Bdissimilarity","Bsecond_moment","Bmean")
    classified<-terra::predict(imgm,classifier_RF,na.rm=T)
    time1 <- proc.time()[3]
    summary$durations[4]<-duration <- as.numeric(time1-time0)
    print(duration)


#------------------------------------------------------------
#  PROCESSING TIME OF RF MODELS INCLUDING BANDS SELECTED BY SFFS
#------------------------------------------------------------
  classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
  bestBands<-summary[summary$name_test=="sel","bands_used"]
  classifier_RF<-readRDS("./Phase_3_RF_classifiers/RF_selected.rds")
  time0 <- proc.time()[3] 
    imgr<-brick("./IMGS_2014/010/SiteJEA010_201404121334.jpg")
    imgr <- terra::aggregate(imgr, df)
    #NB: here compute only the "best Bands selected in SFFS" in our case:
    # rgbvi,gli,vari,ngrdi, R_second_moment, B_contrast, B_entropy and B_second_moment.
    imgr[[4]]<-((imgr[[2]]*imgr[[2]])-(imgr[[1]]*imgr[[3]]))/((imgr[[2]]*imgr[[2]])+(imgr[[1]]*imgr[[3]]))#RGBVI
    imgr[[5]]<-((2*imgr[[2]])-imgr[[1]]-imgr[[3]])/((2*imgr[[2]])+imgr[[1]]+imgr[[3]])#GLI
    imgr[[6]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]]-imgr[[3]])#VARI
    imgr[[6]][imgr[[6]]==Inf]<-max(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))]) #to avoid +inf value in VARI
    imgr[[6]][imgr[[6]]==-Inf]<-min(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))])
    imgr[[7]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]])#NGRDI
    glcm.red <- glcm(imgr[[1]],
                     window = c(ws, ws),
                     na_opt="center",min_x=0,max_x=255,
                     shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                     statistics = c("second_moment"))
    glcm.blue <- glcm(imgr[[3]],
                      window = c(ws, ws),
                      na_opt="center",min_x=0,max_x=255,
                      shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                      statistics = c("contrast",        "entropy",
                                     "second_moment"))
    imgm<-mask(stack(imgr,glcm.red,glcm.blue),roi.data[[1]]$polygons)
    names(imgm)<-c("R","G","B","rgbvi","gli","vari","ngrdi","Rsecond_moment", 
                   "Bcontrast","Bentropy","Bsecond_moment")
    classified<-terra::predict(imgm,classifier_RF,na.rm=T)
    time1 <- proc.time()[3]
    summary$durations[5]<- duration <- as.numeric(time1-time0)
    print(duration)

# write the durations dataframe as a csv file
write.csv(summary,file="./Phase_3_RF_classifiers_accuracies.csv")
```
