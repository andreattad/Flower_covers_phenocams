# 4.1 Flower covers extraction
Once the final RF classifier had been trained, the percentage of pixels in each class was computed for each image. For this, images were selected from the series (see section 2.3.1), for each image the selected features were computed, and percentages of each class within each image was calculated using the RF classifier developed in subsection 2.3.3. 

```r

library(crfsuite)
library(randomForest)
library(glcm)
library(jpeg)
library(raster)
library(terra)
library(lattice)
library(grid)
library(gridExtra)
rm(list=ls())

###-------------------------------------------------------------
###          1)define the feature extraction parameter and import the RF classifier
###------------------------------------------------------------
rf <- readRDS("your/folder/path/Phase_3_RF_classifiers/RF_selected.rds")
ws=11
df=4

bestBands<-names(rf[["forest"]][["xlevels"]])
bestBands

###-------------------------------------------------------------
###          2)define the list of plots with relative block class
###------------------------------------------------------------


plotlist<-c(
  c("001", "002", "003", "005", "007", "008", "009", "010", "011", "013", "016", "017", "018", "020", "021", "024", "025", "026", "027", "028", "030", "033", "035", "037", "039", "040", "042", "043", "044" ,"045", "046"),
  c("048", "049", "051", "053", "054", "056", "057", "058", "059", "060", "062", "064", "065", "067", "070", "071", "073", "074", "075", "077", "079", "080", "081", "082", "083", "084", "085", "088", "090", "091", "092"),
  c("093", "094", "095", "097", "099", "100", "102", "103", "105", "108", "109", "110", "111", "113", "114", "115", "116", "119", "120", "121", "125", "128", "129", "130", "131", "133", "135", "136", "137", "138")
)
blocklist<- c(rep("A",31), rep("B",31) ,rep("C",30))
blockplotlist<-paste(blocklist,plotlist,sep="")
blockplotlist<-blockplotlist[! blockplotlist %in% c("A020",'A033',"A027")] #REMOVED plots in which image acquisition failed.

###-------------------------------------------------------------
###          3) extract a flower cover time series for each plot
###------------------------------------------------------------
#expected processing time= 9 sec for each image, around 137 images per plot, 89 plots
#for each plot->25 minutes per plot.
bp<-blockplotlist[1]
for (bp in blockplotlist){
  time0 <- proc.time()[3];  time0.all <- proc.time()[3]
  block<- substr(bp,1,1)
  plot<- substr(bp,2,4)
    setwd(paste0("your/folder/path/IMGS_2014/",plot))
    tab<-read.csv(paste0("your/folder/path/Phase_1_2014_filtered_indices_BRIAV_BRISD/rawdatafilt",plot,".csv"))
    tab$date<-as.character(as.POSIXct(tab$date,tryFormats = c("%Y-%m-%d %H:%M:%S","%Y-%m-%d %H:%M",
                                                              "%d/%m/%Y %H:%M")))
    tab$imname<-paste("SiteJE",block,plot,"_2014",substr(tab$date,6,7),#yearmonth
                      substr(tab$date,9,10),substr(tab$date,12,13),#dayhour
                      substr(tab$date,15,16),".jpg",sep="")#minute
    datafs<-tab[tab$selSUM==2,] 
    imlist<-datafs$imname
    datafs$Gra_flower<-datafs$Green_vegetation<-datafs$Kna_arv_flower<-datafs$Leu_vul_flower<-datafs$Ran_acr_flower<-datafs$Soil<-NA
    
   i<-1     
   for (i in 1:length(imlist)){
        img<-brick(imlist[i])
          load(paste("your/folder/path/ROISREFS_2014/",block,"/",plot,"/ROI/roi.data.Rdata",sep=""))          
          imgr <- terra::aggregate(img, df)
          #NB: here compute only the "best Bands selected in SFFS" in our case: rgbvi,gli,vari,ngrdi, R_second_moment, B_contrast, B_entropy and B_second_moment.
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
          classified<-terra::predict(imgm,rf,na.rm=T)
          vals0<-summary(as.factor(classified@data@values))
          notNAcount<-sum(subset(vals0,!names(vals0)=="NA's"))
          rf[["classes"]]
          datafs[i, rf[["classes"]][1]]<-if(length(subset(vals0,names(vals0)=="1"))>0){subset(vals0,names(vals0)=="1")/notNAcount}else{NA}
          datafs[i,rf[["classes"]][2]]<-if(length(subset(vals0,names(vals0)=="2"))>0){subset(vals0,names(vals0)=="2")/notNAcount}else{NA}
          datafs[i,rf[["classes"]][3]]<-if(length(subset(vals0,names(vals0)=="3"))>0){subset(vals0,names(vals0)=="3")/notNAcount}else{NA}
          datafs[i,rf[["classes"]][4]]<-if(length(subset(vals0,names(vals0)=="4"))>0){subset(vals0,names(vals0)=="4")/notNAcount}else{NA}
          datafs[i,rf[["classes"]][5]]<-if(length(subset(vals0,names(vals0)=="5"))>0){subset(vals0,names(vals0)=="5")/notNAcount}else{NA}
          datafs[i,rf[["classes"]][6] ]<-if(length(subset(vals0,names(vals0)=="6"))>0){subset(vals0,names(vals0)=="6")/notNAcount}else{NA}
          time1 <- proc.time()[3]
          duration <- time1-time0
          print(duration)
      }
        time1.all <- proc.time()[3];    duration.all <- time1.all-time0.all
        write.csv(datafs,file = paste("your/folder/path/Phase_4_FCTS/dataflowers",plot,".csv",sep=""))
    }



```
