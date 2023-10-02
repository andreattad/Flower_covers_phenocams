# 2.1 Feature values computation
To increase the spectral separability of pixels between different classes, we computed RGB-based features: vegetation indices and texture metrics, described in detail in Table 1. We selected four vegetation indices well established in colour analysis literature (Lussem et al, 2018; Zhao, 2021). For pixels with a specific shade of purple colour, the calculation of the Visible Atmospherically Resistant Index (VARI) resulted in infinite values (for definition, see Table 1). Since only finite values can be used for classifier development, infinite VARI values were replaced with the highest finite value sampled (or lowest in case of negative infinite values), which occurred in less than 0.1% of the labelled pixels. The image textures were derived from co-occurrence matrices for each colour band, since we expected that the flower colours differed from the background (green vegetation or soil) surfaces (Guru et al., 2010), using the “glcm” package in R Studio (Zvoleff, 2020; Haralick et al, 1973). Homogeneity, Contrast, Dissimilarity, Entropy, Second moment, Mean, and Variance were computed in four directions (0°, 45°, 90° and 135°) and then averaged to one rotation-invariant texture as commonly used in texture analysis (e.g., Guru et al., 2010). For the computation of texture metrics, we needed to define the size of the window used for co-occurrence matrices. Moreover, downscaling the images to a lower resolution before feature extraction can give the best detection accuracy while also vastly increasing processing speed compared to higher resolution images (Mann et al., 2022).

We tested the influence of the size of the window used for co-occurrence matrices computation (WS) and downscaling factor (DF) on classification accuracy and processing time. The tested values for WS were 3, 5, 7, 11, 19, 27, and 43 pixels (which is 3, 3+2, 3+4, 3+8, 3+16, 3+24, 3+40) and for DF, we used 1, 2, 4, and 8 pixel. Since flower dimensions were generally much smaller than one tenth of the image height, DF-WS combinations resulting in a window larger than one tenth of the image height were discarded (DF*WS<1040/10). The resulting number of DF-WS combinations was 23.

Here, we compute the 28 feature values for each of the downscaling factor - window size combinations. Prediction capability of random forest models developed based on the extracted features will be compared in phase 2.2
```r
library(dplyr)
library(randomForest)
library(crfsuite)
library(raster)
library(terra)
library(glcm)
setwd(maindir<-"path/to/your/working/directory/")

# 1 Load the labelled pixels dataset
  labs<-read.csv(row.names=1,"./Phase_1_labelled.csv",stringsAsFactors = T)
  imgs<-unique(labs$im)

# 2 Define df-ws combinations and extract the image features in the labelled pixel 
  dfs<-c("08","04","02","01")
  wss<-c("43","27","19","11","07","05","03")
  combinations<-data.frame(dfs=as.numeric(sort(rep(dfs,length(wss)))),wss=as.numeric(rep(wss,length(dfs))))
  combinations$product<-combinations$dfs*combinations$wss
  combinations<-combinations[combinations$product<104,]
  
  for(f in 1:nrow(combinations)){
    df<-combinations$dfs[f];df
    ws<-combinations$wss[f];ws
    
      for (imgname in imgs){

        # Load the image and aggregate it
          img<-brick(as.character(imgname))
          img<-crop(img,extent(0,1280,16,1040)) #16 rows at the bottom of the pictures contain the timestamps.
          imgr <- terra::aggregate(img,fact=df)

        # Load the points and define the minimum extent required to compute points features
          labsi<-labs[labs$im==imgname,]
          pt<- SpatialPoints(coords =labsi[,c("x","y")])
          cellnum<-raster::extract(imgr,pt,sp=T,cellnumbers=T)@data[["cells"]]
          rown<-rowColFromCell(imgr, cellnum)[,1]
          coln<-rowColFromCell(imgr, cellnum)[,2]
          xmin<-((min(coln)-((ws-1)/2))-2)*df
          xmax<-((max(coln)+((ws-1)/2))+1)*df
          ymin<-extent(imgr)[4]-((max(rown)+((ws-1)/2))+1)*df
          ymax<-extent(imgr)[4]-((min(rown)-((ws-1)/2))-2)*df
          new.extent<-extent(xmin,xmax,ymin,ymax)
          cropped<-crop(imgr,new.extent)

        # Compute the image features
          cropped[[4]]<-((cropped[[2]]*cropped[[2]])-(cropped[[1]]*cropped[[3]]))/
                        ((cropped[[2]]*cropped[[2]])+(cropped[[1]]*cropped[[3]]))#RGBVI
          cropped[[5]]<-((2*cropped[[2]])-cropped[[1]]-cropped[[3]])/
                        ((2*cropped[[2]])+cropped[[1]]+cropped[[3]])#GLI
          cropped[[6]]<-(cropped[[2]]-cropped[[1]])/(cropped[[2]]+cropped[[1]]-cropped[[3]])#VARI
          cropped[[7]]<-(cropped[[2]]-cropped[[1]])/(cropped[[2]]+cropped[[1]])#NGRDI
          glcm.red <- glcm(cropped[[1]],
                           window = c(ws, ws),
                           na_opt="center",min_x=0,max_x=255,
                           shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                           statistics = c("variance","homogeneity","contrast","entropy",
                                          "dissimilarity", "second_moment","mean"))
          glcm.green <- glcm(cropped[[2]],
                             window = c(ws, ws),
                             na_opt="center",min_x=0,max_x=255,
                             shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                             statistics = c("variance","homogeneity","contrast","entropy",
                                            "dissimilarity", "second_moment","mean"))
          glcm.blue <- glcm(cropped[[3]],
                            window = c(ws, ws),
                            na_opt="center",min_x=0,max_x=255,
                            shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                            statistics = c("variance","homogeneity","contrast","entropy",
                                           "dissimilarity", "second_moment","mean"))

        # Stack the features, extract their values at labelled pixels
          imgrs<-stack(cropped,glcm.red,glcm.green,glcm.blue)
          ptsa<-raster::extract(imgrs,pt,sp=T)
          ptsa@data$x<-ptsa@coords[,1]
          ptsa@data$y<-ptsa@coords[,2]
          ptsa@data$type<-labsi$type
          ptsa@data$imname<-imgname
          ptsa_data<-ptsa@data
          colnames(ptsa_data)<-c("R","G","B","rgbvi","gli","vari","ngrdi",
               "Rvariance","Rhomogeneity","Rcontrast","Rentropy","Rdissimilarity","Rsecond_moment","Rmean",
               "Gvariance","Ghomogeneity","Gcontrast","Gentropy","Gdissimilarity","Gsecond_moment","Gmean",
               "Bvariance","Bhomogeneity","Bcontrast","Bentropy","Bdissimilarity","Bsecond_moment","Bmean",
               "x","y","type","imname")
   # Stack all extracted features values in a file for each df-ws combination and save it as a csv file
        if(i==1){totsampled<-ptsa_data}else{totsampled<-rbind(totsampled,ptsa_data)}
      }
    write.csv(totsampled,file = paste0("./Phase_2_df_ws_ext_feat/Phase2_df",
                                       as.character(df),"_ws",as.character(ws),".csv"))
  }
```

# 2.2 Model comparison
Here, we tested the influence of window size and downscaling factor on classification accuracy.
We randomly assigned pixels from 70% of the images for training, and 30% for validation for each DF-WS combination, then we trained a random forest (RF) classifier using the “randomForest” function of the “randomForest” package (Liaw & Wiener, 2002). Accuracies were reported and saves as .csv file. 
The metric to calculate the accuracy of the RF classifiers was the mean F1 score of the five classes. The F1 score is derived from precision and recall metrics. The precision is intuitively the ability of the classifier not to label a sampled pixel as positive when it is negative, whereas the recall is the ability of the classifier to find all the positive sampled pixels. Metrics equations are reported in the main text. 
```r
library(randomForest)
library(crfsuite)
library(readr)
library(dplyr)
library(raster)
rm(list=ls())

setwd(maindir<-"path/to/your/working/directory/")
files<-list.files(path="./Phase_2_df_ws_ext_feat/", pattern="Phase2")
res<-data.frame(f1_mean=c(rep(NA,length(files))),
                df=c(rep(NA,length(files))),
                ws=c(rep(NA,length(files))))

for (i in 1:length(files) ){
  print(files[i])

  # Load the dataset and replace VARI infinite values
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
    data$tv<-NULL #We create a new column that labels each observation as either "training" or "validation"
    train<-data[is.element(data$imname,imgst),c("tv")]<-"t"
    valid<-data[!is.element(data$imname,imgst),c("tv")]<-"v"
    train<-data[data$tv=="t",]
    valid<-data[data$tv=="v",]
    summary(train$type); summary(valid$type)

  # Build the classifier
    set.seed(1509)  # Setting seed
    classifier_RF = randomForest(x = train[,c(1:28)],
                               y = as.factor(train$type),
                               ntree = 500,proximity=T)
    pred_rf = predict(classifier_RF, newdata = valid)

  # Calculate its accuracy and report it in the results dataframe
    print(crf_evaluation(pred_rf, valid$type))
    print(table(pred_rf, valid$type))
    res$f1_mean[i]<-crf_evaluation(pred_rf, valid$type)[["overall"]][["f1_mean"]]
    res$df[i]<-parse_number(substr(files[i],10,11))
    res$ws[i]<-parse_number(substr(files[i],14,15))
}
# Write the results dataframe as a csv file
write.csv(res,file="./Phase_2_dfws_accuracy.csv")
```
