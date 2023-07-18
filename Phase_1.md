# 1.1 Brightness and contrast extraction

Light conditions heavily affect pixel colours: images with high brightness were usually foggy, and images with high contrast were usually acquired in direct sunlight conditions. We calculated brightness and contrast for all images using the “extractVIs” function of the R package “Phenopix” (Filippa et al., 2016) and tested which brightness and contrast combinations allow the selection of images acquired in homogeneous light conditions.
Directory structure in the "your/path/folder/ROISREFS_2014" is compatible with the one required in phenopix, a very commonly used package for phenocam images processing. The extraction of each time-series required around 20 minutes.

In this script, our aim is to extract brightness average and contrast for each image.

```r
library(phenopix)

plotIDs<-c("001", "002", "003", "005", "007", "008", "009", "010", "011", "013", "016", "017", "018",  "021", "024", "025", "026", "028", "030", "035", "037", "039", "040", "042", "043", "044" ,"045", "046","048", "049", "051", "053", "054", "056", "057", "058", "059", "060", "062", "064", "065", "067", "070", "071", "073", "074", "075", "077", "079", "080", "081", "082", "083", "084", "085", "088", "090", "091", "092","093", "094", "095", "097", "099", "100", "102", "103", "105", "108", "109", "110", "111", "113", "114", "115", "116", "119", "120", "121", "125", "128", "129", "130", "131", "133", "135", "136", "137", "138")

#######------ Extract vegetation indices -------####
time0.all <- proc.time()[3]
for (plot in plotIDs) {
  time0 <- proc.time()[3]
  path.JPG <- paste("your/folder/path/IMGS_2014/", plot"/", sep="")
  path.ROI <- paste("your/folder/path/ROISREFS_2014/", plot, "/ROI/", sep="")
  path.VI <- paste("your/folder/path/ROISREFS_2014/", plot, "/VI/", sep="")
  setwd(path.JPG);
  extracted<-extractVIs(path.JPG, path.ROI, vi.path=path.VI,plot=F, spatial=F, date.code="yyyymmddHHMM")
  duration <- proc.time()[3]-time0;
  print(duration)
}

# Convert from Rdata format to csv format and save all them in the same folder
for (plot in plotIDs) {
  load(paste0("your/folder/path/ROISREFS_2014/", plot, "/VI/VI.data.Rdata"))
  path.csv <- "your/folder/path/Phase_1_2014_raw_indices/"
  write.csv(VI.data$roi1,file=paste0(path.csv,"VI_raw2014_", plot,".csv"))
}
```
# 1.2 Images selection based on brightness and contrast
Images with uniform light conditions were retrieved by selecting brightness and contrast between the 10th and the 40th percentile within a 3-day window. The selection of the best images within this 3-day window avoided including images taken on days with sub-optimal observations (e.g., all foggy/high contrast images).

```r
#######------DEFINE SETTINGS -------####
start=113           #Start of the season of interest (doy)
end=150             #End of the season of interest (doy)
perc<-c(0.1,0.4)    #Images with uniform light conditions were retrieved by selecting brightness and contrast between the 10th and the 40th percentile in a 3-days window

for (plot in 1:length(plotIDs)){
   dataraw<- read.csv(paste0("your/folder/path/Phase_1_2014_raw_indices/"VI_raw2014_", plot,".csv"),encoding="UTF-8",row.names=1 )
   doys<-unique(dataraw$doy)
   doys<-doys[doys>start&doys<end]

   #######------FIND IMAGES MATCHING BRIGHTNESS CRITERIA-------####
   for (i in 1:length(doys)){
        tab3daysbri<-dataraw[dataraw$doy>doys[i]-2&dataraw$doy<doys[i]+2,]
        tab3daysbri$percentiles<-ecdf(tab3daysbri$bri.av)(tab3daysbri$bri.av)
        dates_selected_bri_i<-tab3daysbri[tab3daysbri$percentiles>perc[1]&tab3daysbri$percentiles<perc[2],])
        if(i==1){dates_selected_bri_tot<-dates_selected_bri_i}else{
                 dates_selected_bri_tot<-rbind(dates_selected_bri_tot,dates_selected_bri_i)
        }
    # eliminate duplicated (some images may be selected more than once)
    dates_selected_bri_tot<-dates_selected_bri_tot[!duplicated(dates_selected_bri_tot$date), ]
    
    # add column selBRI in which images matching the criteria of selection for brightness are TRUE 
    match_bri<-is.element(dataraw$date, dates_selectedtot$date)
    dataraw$selBRI<-FALSE
    dataraw$selBRI[match_bri]<-TRUE

    #######------FIND IMAGES MATCHING CONTRAST CRITERIA-------####
    for (i in 1:length(doys)){
       tab3dayscon<-dataraw[dataraw$doy>doys[i]-2&dataraw$doy<doys[i]+2,]
       tab3dayscon$percentiles<-ecdf(tab3dayscon$bri.sd)(tab3dayscon$bri.sd)
       dates_selected_con_i<-tab3dayscon[tab3dayscon$percentiles>perc[1]&tab3dayscon$percentiles<perc[2],])
        if(i==1){dates_selected_con_tot<-dates_selected_con_i}else{
                 dates_selected_con_tot<-rbind(dates_selected_con_tot,dates_selected_con_i)
        }
    # eliminate duplicated (some images may be selected more than once)
    dates_selected_con_tot<-dates_selected_con_tot[!duplicated(dates_selected_con_tot$date), ]
    
    # add column selBRI in which images matching the criteria of selection for brightness are TRUE 
    match_con<-is.element(dataraw$date, dates_selected_con_tot$date)
    dataraw$selCON<-FALSE
    dataraw$selCON[match_con]<-TRUE


    #######------FIND IMAGES MATCHING BOTH CRITERIA-------####
    dataraw$selBRICON<-FALSE
    dataraw$selBRICON[dataraw$selBRI==T&dataraw$selCON==T]<-T
     
    #######------SAVE THE RESULTS-------####
    block<-if (as.numeric(plot)<48){"A"}else if (as.numeric(plot)<93){"B"}else{"C"}
    dataraw$imname<-paste0("SiteJE",block,plot,"_2014",#year
                                    substr(dataraw$date,6,7),#month
                                    substr(dataraw$date,9,10),#day
                                    substr(dataraw$date,12,13),#hour
                                    substr(dataraw$date,15,16),".jpg")#minute
     dataraw<-dataraw[dataraw$doy>start&dataraw$doy<end,]
     setwd("your/folder/path/Phase_1_2014_filtered_indices_BRI_CON/")
     write.csv(dataraw,file = paste("rawdatafilt",plot,".csv",sep=""))
}

```

# 1.3 Image labelling
Image labelling phase includes: i) preparation of the list of images to be labelled, ii) label the images, iii) remove labelled pixels which are too close from one another. The labelled dataset is provided in the ETH repository (see README file for details) for replicability.

## Image list preparation
To develop a labelled dataset, 300 images were randomly selected (60 images in the period between Apr, 24 and May, 5; 60 images between May, 6 and May, 18; 120 images between May, 19 and May, 29). 
```r
library(stringr)
####---------------1)PREPARE THE LIST OF IMAGES TO BE LABELLED----------------###
setwd("your/folder/path/Phase_1_2014_filtered_indices_BRICON/")
list<-list.files(pattern="filt")

#### CREATE A LIST OF ALL THE FILTERED IMAGES ###
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

#### SAMPLE 60,60 and 180 images in the first, second and third periods, rispectively.###
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
          paste0("C:/Users/david/Documents/dottorato 2020/estero/JENA_2014_MAT_RES/Phase_1_check_sel_imgs/",
                 sapply(longunsorted, function(x) str_extract(x, "(SiteJE\\w+)")))
          )
write.csv(longunsorted,file = "your/folder/path/Phase_1_imgslist300.csv")
```
## Image patches  labelling
For each image, a 200 pixels × 200 pixels image patch was randomly selected and plotted in RGB colours using the “plotRGB()” function of the “raster” package (Hijmans, 2022). Around 30 pixels per image were labelled by clicking on the image to retrieve the x and y coordinates using the “locator()” function of the “graphics” package and assigning to each pixel the class to which it belongs (see subsection 2.1). 

```r

library(raster)
library(rgdal)
rm(list=ls());

# load the image list and create seeds list that you will use to select a random 200X200 detail in the image
imgslist<-read.csv(file = "your/folder/path/Phase_1_imgslist300.csv",
                      row.names=1)[,1]
set.seed(1002); seedsx<-sample(1:2^15,300)
set.seed(1056); seedsy<-sample(1:2^15,300)
buffer<-100

####  LABEL PIXEL OF CLASS SOIL IN THE PICTURES
####  1) define the classes
classes <- c("Soil", "Green_vegetation", "Kna_arv_flower", "Ran_acr_flower", "Leu_vul_flower", "Gra_flower")
####  2) select one class (e.g, Soil)
class<-classes[1]
####  3) run line 76, an image will appear in the plot pane.
####  4) click on some pixels of the class selected in point 2) (e.g., Soil)(if there are some), then press esc. A dataframe with coordinates, image name and label will be saved as a ".csv" file.
####  5) press Ctrl+Alt+T to run the section again. A new image will appear in the plot pane. Repeat point 4) and 5) until all the images have been labelled
####  6) when you will have put some labels on all images, go to point 2) and select the second class. Then go on with point 3),4),5),6) until all classes will be labelled.  
####  NB: LOCATOR WORKS PROPERLY ONLY WHEN ZOOM IN GLOBAL OPTIONS IS SET TO 100 %
   

####--------------------- RUN HERE
{
   if(!exists("i")) {i <- 1} else {i<-i+1};  print(paste0("i=",i,"class=",class))
    img<-zoiCx<-zoiCy<-NULL
   img<-brick(imgslist[i])
   set.seed(seedsx[i]);  zoiCx<-sample(200:1080,1);
   set.seed(seedsy[i]);  zoiCy<-sample(200:840,1)
   imgc <- crop(img, extent(zoiCx-buffer,zoiCx+buffer,zoiCy-buffer,zoiCy+buffer))
   par(xaxt = "n",yaxt = "n",mar=c(0,0,0,0),oma=c(0,0,0,0))
   plotRGB(imgc,axes = F)
   labelled<-NULL
   labelled<-data.frame(locator(type="p",pch=16,col="red"))
   if(nrow(labelled)>0) { 
      labelled$type<-class; labelled$im<-imgslist[i];labelled$zoiCx<-zoiCx;labelled$zoiCy<-zoiCy
   write.csv(labelled,file=paste0("your/folder/path/Phase_1_lab_xy/",
                                "lab_xy_pic_",sprintf("%03d", i),"class_",class,".csv"))}
    }

```

## Close pixels removal
To prevent duplicated pixels after downscaling the images, any labelled pixels that were within a distance of 8 pixels from one another were eliminated from the dataset. The labelling phase resulted in a table where the class and pixel coordinates were stored.
```r
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
write.csv(ProxFiltered_tot,"your/folder/path/Phase_1_labelled.csv")
```


