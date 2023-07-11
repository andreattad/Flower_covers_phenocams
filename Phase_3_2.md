# 3.2 Processing time calculation






```r
library(randomForest)
library(glcm)
library(raster)
library(terra)
library(crfsuite)
rm(list=ls())

summary<-read.csv("your/folder/path/Phase_3_RF_classifiers_accuracies.csv",row.names=1)
summary$durations<-c(NA,NA,NA,NA,NA)

df<-4
ws<-11
load(paste("your/folder/path/ROISREFS_2014/A/010/ROI/roi.data.Rdata",sep=""))


#------------------------------------------------------------
#  1 PROCESSING TIME OF RF MODELS INCLUDING ALL BANDS
#------------------------------------------------------------
classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
classifier_RF<-readRDS("your/folder/path/Phase_3_RF_classifiers/RF_all.rds")
  time0 <- proc.time()[3]
      imgr<-brick("your/folder/path/IMGS_2014/010/SiteJEA010_201404121334.jpg")
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
      names(imgm)<-c("R","G","B","rgbvi","gli","vari","ngrdi",
                             "Rvariance","Rhomogeneity","Rcontrast","Rentropy","Rdissimilarity","Rsecond_moment","Rmean",
                             "Gvariance","Ghomogeneity","Gcontrast","Gentropy","Gdissimilarity","Gsecond_moment","Gmean",
                             "Bvariance","Bhomogeneity","Bcontrast","Bentropy","Bdissimilarity","Bsecond_moment","Bmean")
      # plot(is.na(imgm))
      classified<-terra::predict(imgm,classifier_RF,na.rm=T)
      time1 <- proc.time()[3]
      summary$durations[1]<-duration <- as.numeric(time1-time0)
      print(duration)
      row.names(classifier_RF[["importance"]])
      
      #------------------------------------------------------------
      #  2 PROCESSING TIME OF RF MODELS INCLUDING RGB BANDS
      #------------------------------------------------------------
      classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
      
      classifier_RF<-readRDS("your/folder/path/Phase_3_RF_classifiers/RF_rgb.rds")
      time0 <- proc.time()[3]
      imgr<-brick("your/folder/path/IMGS_2014/010/SiteJEA010_201404121334.jpg")
      imgr <- terra::aggregate(imgr, df)
      imgm<-mask(imgr,roi.data[[1]]$polygons)
      names(imgm)<-c("R","G","B")
      classified<-terra::predict(imgm,classifier_RF,na.rm=T)
      time1 <- proc.time()[3]
      summary$durations[2]<-duration <-  as.numeric(time1-time0)
      print(duration)
      row.names(classifier_RF[["importance"]])
      
 
      
      
      #------------------------------------------------------------
      #  3 PROCESSING TIME OF RF MODELS INCLUDING RGB+VI
      #------------------------------------------------------------
      classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
      
      classifier_RF<-readRDS("your/folder/path/Phase_3_RF_classifiers/RF_rgb_vi.rds")
      time0 <- proc.time()[3]
      imgr<-brick("your/folder/path/IMGS_2014/010/SiteJEA010_201404121334.jpg")
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
      row.names(classifier_RF[["importance"]])
      #------------------------------------------------------------
      #  4 PROCESSING TIME OF RF MODELS INCLUDING RGB+TEXTURE METRICS
      #------------------------------------------------------------
      classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL

      classifier_RF<-readRDS("your/folder/path/Phase_3_RF_classifiers/RF_rgb_texm.rds")
      time0 <- proc.time()[3]
      imgr<-brick("your/folder/path/IMGS_2014/010/SiteJEA010_201404121334.jpg")
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
      row.names(classifier_RF[["importance"]])
     
      
      #------------------------------------------------------------
      #  PROCESSING TIME OF RF MODELS INCLUDING BANDS SELECTED BY SFFS
      #------------------------------------------------------------
      classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
      
      bestBands<-summary[summary$name_test=="sel","bands_used"]
      bestBands
      classifier_RF<-readRDS("your/folder/path/Phase_3_RF_classifiers/RF_selected.rds")
      time0 <- proc.time()[3]
      imgr<-brick("your/folder/path/IMGS_2014/010/SiteJEA010_201404121334.jpg")
      imgr <- terra::aggregate(imgr, df)
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
      row.names(classifier_RF[["importance"]])
      write.csv(summary,file="your/folder/path/Phase_3_RF_classifiers_accuracies.csv")
      dev.off()
      par(mfrow=c(2,1))
      text(x=barplot(summary$f1, names.arg = summary$name_test),y=0.2, label=round(summary$f1,3),adj=c(0.5,1))
      text(x=barplot(summary$durations, names.arg = summary$name_test),y=3, label=round(summary$durations),adj=c(0.5,1))
      
      ```
