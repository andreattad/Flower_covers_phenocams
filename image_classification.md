# IMAGE CLASSIFICATION
Here, we will classify and save the images of one plot. 
![Example of RGB and classified image from Figure 5 in the manuscript](https://drive.google.com/uc?id=1NVcvDAzGoqVIJ4gtlL2xgXSbHAvd22ZY)  
*Example of RGB and classified image from Figure 5 in the manuscript*

```r
library(randomForest)
library(raster)
library(glcm)
library(terra)
library(lattice)
library(grid)
library(gridExtra)
library(stringr)

#---------upload the RF classifier
rm(list=ls())
setwd("your/folder/path/")
rf <- readRDS("your/folder/path/Phase_3_RF_classifiers/RF_selected.rds")
#----------define the plot you want to classify
plot<-"056"
#----------read the CSV containing the selected images IDs based on brightness and contrast
csvlist<-list.files("Phase_1_2014_filtered_indices_BRIAV_BRISD",full.names=T)
csvplot<-csvlist[grepl(plot,csvlist)]
imlist<-read.csv(csvplot)
imlist<-imlist[imlist$selSUM==2,"imname"]
fullpathimlist<-paste0("IMGS_2014/",plot,"/",imlist)

#---------create a folder where the classified images will be stored
dir<-paste0("your/folder/path/Phase_4_Classified_images_various/",plot)
dir.create(dir)

#----------define the ws and df, upload the images, compute the features, classify the image and save a false color plot of the classified image
ws=11
df=4
i<-3
for (i in 1:length(fullpathimlist)){
   imgn<-fullpathimlist[i]
   time0 <- proc.time()[3];  time0.all <- proc.time()[3]
   
   img<-brick(imgn)
   imgr <- terra::aggregate(img, df)
 #NB: here compute only the "best Bands selected in SFFS" in our case: rgbvi,gli,vari,ngrdi, R_second_moment, B_contrast, B_entropy and B_second_moment
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
   imgm<-stack(imgr,glcm.red,glcm.blue)
   names(imgm)<-c("R","G","B","rgbvi","gli","vari","ngrdi","Rsecond_moment", 
                  "Bcontrast","Bentropy","Bsecond_moment")
   classified<-terra::predict(imgm,rf,na.rm=T) 
   newname<-paste0(dir,"/",imlist[i])
   
   tiff(filename = paste0(newname,"classified.tiff"),width=1280,height=1040)
   rf[["classes"]]#this is the order of values
   pcla<-rasterVis::levelplot(classified,  margin=F,colorkey=F,scales=list(draw=F),
                              at = seq(1,6), col.regions=c("grey20","seagreen","purple","cyan","gold","brown"))
   
   plot(pcla)
   dev.off()
   time1 <- proc.time()[3]
   duration <- time1-time0
   print(duration)
}

#----------------Copy the original images in the same subfolder to compare classified and RGB
for (i in 1:length(imlist)) {
   file.copy(from=fullpathimlist[i], 
             to=  paste0(dir,"/",imlist[i]))
}
```
