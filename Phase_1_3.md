# Image labelling
Image labelling phase includes: i) prepare the list of images to be labelled, ii) label the images, iii) remove labelled pixels which are too close from one another 

## Image list preparation
To develop a labelled dataset, 300 images were randomly selected (60 images in the period between Apr, 24 and May, 5; 60 images between May, 6 and May, 18; 120 images between May, 19 and May, 29). 
```
####---------------1)PREPARE THE LIST OF IMAGES TO BE SAMPLED----------------###
setwd("your/folder/path/Phase_1_2014_filtered_indices_BRIAV_BRISD/")
list<-list.files(pattern="filt")
blocks<-c(rep("A",28),rep("B",31),rep("C",30))#should be 31,31,30 but I removed plot 20,27,33, all in block A
####CREATE LIST OF ALL THE FILTERED IMAGES###
for (i in 1:length(list)){
  plot<-substr(list[i],12,14)
  block<-blocks[i]
  tab<-read.csv(list[i],row.names=1)
  datafs<-tab[tab$selSUM==2,]
     if(plot=="001"){  imlist<-datafs$imname} else{
      imlist<-c(imlist,datafs$imname) }
     }

imlistpath<-paste0("your/folder/path/IMGS_2014/",substr(imlist,8,10),"/",imlist)
imlistpath[1]


set.seed(1536)
listperiod1<-sample( grep("_20140424|_20140425|_20140426|_20140427|_20140428|_20140429|_20140430|_20140501
        |_20140502|_20140503|_20140504|_20140505", imlistpath,value=T), 60);set.seed(1536)
listperiod2<-  sample(grep("_20140506|_20140507|_20140508|_20140509|_20140510|_20140511|_20140512|_20140513
      |_20140514|_20140515|_20140516|_20140517|_20140518", imlistpath,value=T), 60);set.seed(1536)
listperiod3<-sample(grep("_20140519|_20140520  |_20140521|_20140522|_20140523|_20140524|_20140525|
                      _20140526|_20140527|_20140528|_20140529", imlistpath,value=T), 180)

longsortedlist<-c(listperiod1,listperiod2,listperiod3)
longsortedlist[101]
set.seed(1536)
longunsorted<-sample(longsortedlist,c(60+60+180))
longunsorted[length(longunsorted)]

file.copy(longunsorted,
          paste0("your/folder/path/Phase_1_check_sel_imgs/",substr(longunsorted,80,120))
          )
write.csv(longunsorted,file = "your/folder/path/Phase_1_imgslist300.csv")
```
## Images patches  labelling
For each image, a 200 pixels × 200 pixels image patch was randomly selected and plotted in RGB colours using the “plotRGB()” function of the “raster” package (Hijmans, 2022). Around 30 pixels per image were labelled by clicking on the image to retrieve the x and y coordinates using the “locator()” function of the “graphics” package and assigning to each pixel the class to which it belongs (see subsection 2.1). 

```

library(raster)
library(rgdal)
rm(list=ls());

#loadimagelist and create seeds list that you will use to select a random 200X200 detail in the image
imgslist<-read.csv(file = "your/folder/path/Phase_1_imgslist300.csv",
                      row.names=1)[,1]
imgslist
set.seed(1002); seedsx<-sample(1:2^15,300)
set.seed(1056); seedsy<-sample(1:2^15,300)
buffer<-100
dev.off()


####-----now define the classes and create a directory for each class--------------------###
classes <- c("Soil", "Green_vegetation", "Kna_arv_flower", "Ran_acr_flower", "Leu_vul_flower", "Gra_flower")


####  LABEL PIXEL OF CLASS SOIL IN THE PICTURES
####  1) run line 99
####  2) click the soil pixel in the window
####  3) press Ctrl+Alt+T to run the section again, until you labelled all the pictures
class<-classes[6]

####---------------------
{
   if(!exists("i")) {i <- 1} else {i<-i+1};  print(paste0("i=",i,"class=",class))
    img<-zoiCx<-zoiCy<-NULL
   img<-brick(imgslist[i])
   #define the center of the detail region you will sample based on random sample
   #we don't want to be label image edges to avoid NA values of features when downscaling and computing texture metrics
   set.seed(seedsx[i]);  zoiCx<-sample(200:1080,1);
   set.seed(seedsy[i]);  zoiCy<-sample(200:840,1)
   imgc <- crop(img, extent(zoiCx-buffer,zoiCx+buffer,zoiCy-buffer,zoiCy+buffer))
   par(xaxt = "n",yaxt = "n",mar=c(0,0,0,0),oma=c(0,0,0,0))
   plotRGB(imgc,axes = F)
   labelled<-NULL
   #NB1: LOCATOR WORKS PROPERLY ONLY WHEN ZOOM IN GLOBAL OPTIONS IS SET TO 100
   labelled<-data.frame(locator(type="p",pch=16,col="red"))
   if(nrow(labelled)>0) { 
      labelled$type<-class; labelled$im<-imgslist[i];labelled$zoiCx<-zoiCx;labelled$zoiCy<-zoiCy
   write.csv(labelled,file=paste0("your/folder/path/Phase_1_lab_xy/",
                                "lab_xy_pic_",sprintf("%03d", i),"class_",class,".csv"))}
    }

```

## Remove too close pixels
To prevent duplicated pixels after downscaling the images (refer to the subsequent section for downscaling details), any labelled pixels that were within a distance of 8 pixels from one another were eliminated from the dataset. The labelling phase resulted in a table where the class and pixel coordinates were stored.
```
library(raster)
library(sp)
# combine all data frames into a single data frame
file_paths<-list.files(path="your/folder/path/Phase_1_lab_xy/",
                       full.names=T,pattern=".csv")
combined_df <- do.call(rbind, lapply(file_paths, read.csv,row.names=1)
summary(factor(combined_df$type))
images_ids<-unique(combined_df$im)
minimum_distance<-8

for(i in 1:length(images_ids)){
   df<-combined_df[combined_df$im==images_ids[i],]
   df$progressive<-1:nrow(df)
   ProxFiltered<- df[1,]
   for (j in 2:nrow(df  )){
      pts<- as.matrix(  ProxFiltered[,c("x","y")])
       pt<-  as.numeric(df[j,c("x","y")])
      dists<- spDistsN1(pts, pt, longlat = FALSE)
      exceed<- any(dists<minimum_distance)
      if (exceed==FALSE){
         ProxFiltered<-rbind(  ProxFiltered, df[j,])
      }}

   if(!exists("ProxFiltered_tot")) {ProxFiltered_tot <- ProxFiltered} else {
      ProxFiltered_tot<-rbind(ProxFiltered_tot,ProxFiltered)};  
 }
write.csv(ProxFiltered_tot,"your/folder/path/1_labelled.csv")
```







