# Feature values computation
To increase the spectral separability of pixels between different classes, we computed RGB-based features: vegetation indices and texture metrics, described in detail in Table 1. We selected four vegetation indices well established in colour analysis literature (Lussem et al, 2018; Zhao, 2021). For pixels with a specific shade of purple colour, the calculation of the Visible Atmospherically Resistant Index (VARI) resulted in infinite values (for definition, see Table 1). Since only finite values can be used for classifier development, infinite VARI values were replaced with the highest finite value sampled (or lowest in case of negative infinite values), which occurred in less than 0.1% of the labelled pixels. The image textures were derived from co-occurrence matrices for each colour band, since we expected that the flower colours differed from the background (green vegetation or soil) surfaces (Guru et al., 2010), using the “glcm” package in R Studio (Zvoleff, 2020; Haralick et al, 1973). Homogeneity, contrast, dissimilarity, entropy, second moment, mean, and variance were computed in four directions (0°, 45°, 90° and 135°) and then averaged to one rotation-invariant texture as commonly used in texture analysis (e.g., Guru et al., 2010). For the computation of texture metrics, we needed to define the size of the window used for co-occurrence matrices. Moreover, downscaling the images to a lower resolution before feature extraction can give the best detection accuracy while also vastly increasing processing speed compared to higher resolution images (Mann et al., 2022).

We tested the influence of the size of the window used for co-occurrence matrices computation (WS) and downscaling factor (DF) on classification accuracy and speed. The tested values for WS were 3, 5, 7, 11, 19, 27, and 43 pixels (which is 3, 3+2, 3+4, 3+8, 3+16, 3+24, 3+40) and for DF, we used 1, 2, 4, and 8 pixel. Since flower dimensions were generally much smaller than one tenth of the image height, DF-WS combinations resulting in a window larger than one tenth of the image height were discarded (DF*WS<1040/10). The resulting number of DF-WS combinations was 23.

Here, we compute the 28 feature values for each of the downscaling factor - window size combinations. Prediction capability of random forest models developed based on the extracted features will be compared in phase 2.2
```
library(dplyr)
library(randomForest)
library(crfsuite)
library(raster)
library(terra)
library(glcm)

rm(list=ls())
#------------------------------------------------------------
#  1 LOAD LABEL DATASET
#------------------------------------------------------------
labs<-read.csv(row.names=1,"your/folder/path/Phase_1_labelled.csv",stringsAsFactors = T)
head(labs)
ims<-unique(labs$im)


#-----------------------------------------------------------------------------------
#  2 EXTRACT IMAGE FEATURES VALUES IN THE LABELLED PIXELS USING VARIOUS DF-WS COMBINATIONS
#-----------------------------------------------------------------------------------

#define downscaling levels
dfs<-c("08","04","02","01")
 
#define window size levels
wss<-c("03","05","07","11","19","27","43")

combi0<-data.frame(dfs=as.numeric(sort(rep(dfs,length(wss)))),wss=as.numeric(rep(wss,length(dfs))))
combi0$product<-combi0$dfs*combi0$wss
combi<-combi0[combi0$product<104,]

for(f in 1:nrow(combi)){

  df<-combi$dfs[f];df
  ws<-combi$wss[f];ws
  
    for (i in 1: length(ims)){
      print(i); #  print(ims[i])
      labsi<-ptsa<-cropped<-img<-glcm.blue<-glcm.green<-glcm.red<-pt<-ptsa_data<-rown<-coln<-cellnum<-xmin<-xmax<-ymin<-ymax<-new.extent<-NA
      
          img0<-brick(as.character(ims[i]))
          img<-crop(img0,extent(0,1280,16,1040))
          imgr <- terra::aggregate(img,fact=df)
          labsi<-labs[labs$im==ims[i],]
          pt<- SpatialPoints(coords =labsi[,c("x","y")])
          cellnum<-raster::extract(imgr,pt,sp=T,cellnumbers=T)@data[["cells"]]
          rown<-rowColFromCell(imgr, cellnum)[,1]
          coln<-rowColFromCell(imgr, cellnum)[,2]
          
          xmin<-((min(coln)-((ws-1)/2))-2)*df
          xmax<-((max(coln)+((ws-1)/2))+1)*df
          ymin<-extent(imgr)[4]-((max(rown)+((ws-1)/2))+1)*df#rownumbering from top to bottom but ycoord from bottom to top
          ymax<-extent(imgr)[4]-((min(rown)-((ws-1)/2))-2)*df
          
          new.extent<-extent(xmin,xmax,ymin,ymax)
          
          cropped<-crop(imgr,new.extent)
          cropped[[4]]<-((cropped[[2]]*cropped[[2]])-(cropped[[1]]*cropped[[3]]))/((cropped[[2]]*cropped[[2]])+(cropped[[1]]*cropped[[3]]))#RGBVI
          cropped[[5]]<-((2*cropped[[2]])-cropped[[1]]-cropped[[3]])/((2*cropped[[2]])+cropped[[1]]+cropped[[3]])#GLI
          cropped[[6]]<-(cropped[[2]]-cropped[[1]])/(cropped[[2]]+cropped[[1]]-cropped[[3]])#VARI
          cropped[[7]]<-(cropped[[2]]-cropped[[1]])/(cropped[[2]]+cropped[[1]])#NGRDI

          
          glcm.red <- glcm(cropped[[1]],
                           window = c(ws, ws),
                           na_opt="center",min_x=0,max_x=255,
                           shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                           statistics = c("variance",       "homogeneity","contrast",        "entropy",
                                          "dissimilarity", "second_moment","mean"))
          glcm.green <- glcm(cropped[[2]],
                             window = c(ws, ws),
                             na_opt="center",min_x=0,max_x=255,
                             shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                             statistics = c("variance",       "homogeneity","contrast",        "entropy",
                                            "dissimilarity", "second_moment","mean"))
          glcm.blue <- glcm(cropped[[3]],
                            window = c(ws, ws),
                            na_opt="center",min_x=0,max_x=255,
                            shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                            statistics = c("variance",       "homogeneity","contrast",        "entropy",
                                           "dissimilarity", "second_moment","mean"))
     
               imgr2<-stack(cropped,glcm.red,glcm.green,glcm.blue)
          ptsa<-raster::extract(imgr2,pt,sp=T)
          ptsa@data$x<-ptsa@coords[,1]
          ptsa@data$y<-ptsa@coords[,2]
          ptsa@data$type<-labsi$type
          ptsa@data$imname<-ims[i]
          ptsa_data<-ptsa@data

          colnames(ptsa_data)<-c("R","G","B","rgbvi","gli","vari","ngrdi",
                                     "Rvariance","Rhomogeneity","Rcontrast","Rentropy","Rdissimilarity","Rsecond_moment","Rmean",
                                     "Gvariance","Ghomogeneity","Gcontrast","Gentropy","Gdissimilarity","Gsecond_moment","Gmean",
                                     "Bvariance","Bhomogeneity","Bcontrast","Bentropy","Bdissimilarity","Bsecond_moment","Bmean",
                                     "x","y","type","imname")
          
          if(i==1){totsampled<-ptsa_data
          }else{totsampled<-rbind(totsampled,ptsa_data)}
    }
  write.csv(totsampled,file = paste0("your/folder/path/Phase_2_df_ws_ext_feat/Phase2_df",
                                     as.character(df),"_ws",as.character(ws),".csv"))
}
```





