# **Tutorial for the manuscript "Extracting single species flowering phenology from grassland species mixtures using time-lapse cameras"**
by D.Andreatta, C. Bachofen, M. Dalponte, V. Klaus, N.Buchmann.   
Submitted to Remote Sensing of Environment. 
Contact information: davide.andreatta@phd.unipd.it.

In this tutorial single species flowering phenology time-series and phenological metrics from time-lapse camera images of grasslands captured at the Jena trait-based experiment (Germany) were extracted. 

To ensure that the results of the manuscript can be easily replicated, we decided to base our tutorial on the provided case study and provide code that can be immediately applied to the supplementary materials, including images, labels, and regions of interest, which will be made available on the ETH repository (link to be provided after acceptance). Although the presented workflow is tailored to the specific case study, it can be adapted for use with different case studies with minor modifications. To apply the workflow to new case studies, adjustments should be made to plot IDs, to the period of interest, and to flower class names.

The code is developed in R version 4.3.0 (2023-04-21 ucrt).

<figure>
<img src="https://drive.google.com/uc?id=1z9yc0Tht4S425txHhc42OpQdz4N7g6I6" width="800">
</figure>

_Structure of the proposed workflow (Figure 2 in the manuscript)_

The proposed workflow can be divided in four phases: 
1. [image selection and pixel labelling](Phase_1.md)
     -    1.1 brightness and contrast extraction
     -    1.2 images selection based on brightness and contrast
     -    1.3 image labelling
2. [selection of downscaling factor and window size and feature computation](Phase_2.md)
     -    2.1 feature values extraction
     -    2.2 model comparison
3. [feature selection and final classifier compilation](Phase_3.md)
     -    3.1 feature selection and model comparison
     -    3.2 processing time calculation
4. [FCTS extraction and smoothing and phenological metrics identification.](Phase_4.md)
     -    4.1 FCTS extraction
     -    4.2 FCTS smoothing and display
     -    4.3 Phenological metric identification
  
Moreover, [here](image_classification.md) we provide the code to classify all the images in a plot
![Example of RGB and classified image from Figure 5 in the manuscript](https://drive.google.com/uc?id=1NVcvDAzGoqVIJ4gtlL2xgXSbHAvd22ZY)  
*Example of RGB and classified image*

# To get ready for the tutorial:

## 1. Download the materials:
The materials can be downloaded at the link: *link to be provided after acceptance*.
Materials include: Images, Region of interests and labels.

## 2. Set up the working environment. 

```
setwd("your/folder/path")
dir.create("Phase_1_raw_indices")
dir.create("Phase_1_filtered_indices_BRI_CON")
dir.create("Phase_1_lab_xy")
dir.create("Phase_2_df_ws_ext_feat")
dir.create("Phase_3_RF_classifiers")
dir.create("Phase_4_FCTS")
dir.create("Phase_4_Classified_images")
```

## 3. Install the required R packages

```
packages <- c("crfsuite", "data.table", "dplyr", "forecast", "ggplot2", 
              "ggpubr", "ggtext", "glcm", "grid", "phenopix", "randomForest",
              "raster", "readr", "rgdal", "sp", "stringr", "terra", "varSel", "zoo")
# Load or install missing packages
lapply(packages, function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
})

```


# 1.1 Brightness and contrast extraction

Light conditions heavily affect pixel colours: images with high brightness were usually foggy, and images with high contrast were usually acquired in direct sunlight conditions. We calculated brightness and contrast for all images using the “extractVIs” function of the R package “Phenopix” (Filippa et al., 2016) and tested which brightness and contrast combinations allow the selection of images acquired in homogeneous light conditions.
Directory structure in the "your/path/folder/ROISREFS" is compatible with the one required in phenopix, a very commonly used package for phenocam images processing. The extraction of each time-series required around 20 minutes.

In this script, our aim is to extract brightness average and contrast for each image.

```r
library(phenopix)
setwd("your/folder/path")

plotIDs<-c( "001","002","003","005","007","008","009","010","011","013","016","017","018","021",
"024","025","026","028","030","035","037","039","040","042","043","044","045","046","048","049",
"051","053","054","056","057","058","059","060","062","064","065","067","070","071","073","074",
"075","077","079","080","081","082","083","084","085","088","090","091","092","093","094","095",
"097","099","100","102","103","105","108","109","110","111","113","114","115","116","119","120",
"121","125","128","129","130","131","133","135","136","137","138")

# Extract vegetation indices
time0.all <- proc.time()[3]
for (plot in plotIDs) {
  time0 <- proc.time()[3]
  path.JPG <- paste("./IMGS/", plot"/", sep="")
  path.ROI <- paste("./ROISREFS/", plot, "/ROI/", sep="")
  path.VI <- paste("./ROISREFS/", plot, "/VI/", sep="")
  extracted<-extractVIs(path.JPG, path.ROI, vi.path=path.VI,plot=F, spatial=F, date.code="yyyymmddHHMM")
  duration <- proc.time()[3]-time0;
  print(duration)
}

# Convert from Rdata format to csv format and save all them in the same folder
for (plot in plotIDs) {
  load(paste0("./ROISREFS/", plot, "/VI/VI.data.Rdata"))
  path.csv <- "./Phase_1_raw_indices/"
  write.csv(VI.data$roi1,file=paste0(path.csv,"VI_raw_", plot,".csv"))
}
```
# 1.2 Images selection based on brightness and contrast
Images with uniform light conditions were retrieved by selecting brightness and contrast between the 10th and the 40th percentile within a 3-day window. The selection of the best images within this 3-day window avoided including images taken on days with sub-optimal observations (e.g., all foggy/high contrast images).

```r
setwd("your/folder/path")
# Define settings
start=113           #Start of the season of interest (doy)
end=150             #End of the season of interest (doy)
perc<-c(0.1,0.4)    #Images with uniform light conditions

for (plot in 1:plotIDs){
   # Load the dataframe with brightness and contrast information 
   dataraw<- read.csv(paste0("./Phase_1_raw_indices/VI_raw_", plot,".csv"),
                      encoding="UTF-8",row.names=1 )
   doys<-unique(dataraw$doy)
   doys<-doys[doys>start&doys<end]

   # Find images matching brightness criteria
             for (i in 1:length(doys)){
                  tab3daysbri<-dataraw[dataraw$doy>doys[i]-2&dataraw$doy<doys[i]+2,]
                  tab3daysbri$percentiles<-ecdf(tab3daysbri$bri.av)(tab3daysbri$bri.av)
                  dates_selected_bri_i<-tab3daysbri[tab3daysbri$percentiles>perc[1]&
                                        tab3daysbri$percentiles<perc[2],]
                  if(i==1){dates_selected_bri_tot<-dates_selected_bri_i}else{
                           dates_selected_bri_tot<-rbind(dates_selected_bri_tot,dates_selected_bri_i)
                  }}
              # eliminate duplicated (some images may be selected more than once)
              dates_selected_bri_tot<-dates_selected_bri_tot[!duplicated(dates_selected_bri_tot$date), ]
              
              # add column selBRI in which images matching the criteria of selection for brightness are TRUE 
              match_bri<-is.element(dataraw$date, dates_selected_bri_tot$date)
              dataraw$selBRI<-FALSE
              dataraw$selBRI[match_bri]<-TRUE

    # Find images matching contrast criteria
              for (i in 1:length(doys)){
                 tab3dayscon<-dataraw[dataraw$doy>doys[i]-2&dataraw$doy<doys[i]+2,]
                 tab3dayscon$percentiles<-ecdf(tab3dayscon$bri.sd)(tab3dayscon$bri.sd)
                 dates_selected_con_i<-tab3dayscon[tab3dayscon$percentiles>perc[1]&
                                        tab3dayscon$percentiles<perc[2],]
                  if(i==1){dates_selected_con_tot<-dates_selected_con_i}else{
                           dates_selected_con_tot<-rbind(dates_selected_con_tot,dates_selected_con_i)
                  }}
              # eliminate duplicated (some images may be selected more than once)
              dates_selected_con_tot<-dates_selected_con_tot[!duplicated(dates_selected_con_tot$date), ]
              
              # add column selBRI in which images matching the criteria of selection for brightness are TRUE 
              match_con<-is.element(dataraw$date, dates_selected_con_tot$date)
              dataraw$selCON<-FALSE
              dataraw$selCON[match_con]<-TRUE


    # Find images matching both criteria
              dataraw$selBRICON<-FALSE
              dataraw$selBRICON[dataraw$selBRI==T&dataraw$selCON==T]<-T
     
    # Save the results
              block<-if (as.numeric(plot)<48){"A"}else if (as.numeric(plot)<93){"B"}else{"C"}
              dataraw$imname<-paste0("SiteJE",block,plot,"_2014",#year
                                              substr(dataraw$date,6,7),#month
                                              substr(dataraw$date,9,10),#day
                                              substr(dataraw$date,12,13),#hour
                                              substr(dataraw$date,15,16),".jpg")#minute
               dataraw<-dataraw[dataraw$doy>start&dataraw$doy<end,]
               write.csv(dataraw,file = paste("./Phase_1_filtered_indices_BRI_CON/rawdatafilt",plot,".csv",sep=""))
}

```

# 1.3 Image labelling
Image labelling phase includes: i) image list preparation, ii) image patches labelling, iii) close pixels removal. The labelled dataset is provided in the ETH repository (see README file for details) for replicability.

## Image list preparation
To develop a labelled dataset, 300 images were randomly selected (60 images in the period between Apr, 24 and May, 5; 60 images between May, 6 and May, 18; 120 images between May, 19 and May, 29). 
```r
library(stringr)

setwd("your/folder/path")

# Create a list of all the filtered images
csv_name_list<-list.files(path="./Phase_1_filtered_indices_BRI_CON/",pattern="rawdatafilt")
combined.df <- do.call(rbind , lapply(csv_name_list, read.csv, row.names = 1))
combined.df_sel<-combined.df[combined.df$selBRICON==T,]
combined.df_sel$imlistpath<-paste0(maindir,"IMGS/",
                        substr(combined.df_sel$imname,8,10),"/",combined.df_sel$imname)
write.csv(combined.df_sel,file="./Phase_1_selected_images.csv")

# Identify which images were selected in each of the three periods
listperiod1<-combined.df_sel[as.Date(combined.df_sel$date)>=as.Date("2014-04-24")&
                        as.Date(combined.df_sel$date)<=as.Date("2014-05-05"),"imlistpath"]
listperiod2<-combined.df_sel[as.Date(combined.df_sel$date)>=as.Date("2014-05-06")&
                        as.Date(combined.df_sel$date)<=as.Date("2014-05-18"),"imlistpath"]
listperiod3<-combined.df_sel[as.Date(combined.df_sel$date)>=as.Date("2014-05-19")&
                        as.Date(combined.df_sel$date)<=as.Date("2014-05-29"),"imlistpath"]
# Sample 60,60 and 180 images in the first, second and third periods, respectively
set.seed(1536)
listperiod1shuffled<- sample(listperiod1,60);set.seed(1536)
listperiod2shuffled<- sample(listperiod2,60);set.seed(1536)
listperiod3shuffled<- sample(listperiod3,180);set.seed(1536)
list3periods<-c(listperiod1shuffled,listperiod2shuffled,listperiod3shuffled)
set.seed(1536)
list3periods_shuffled<-sample(list3periods,c(60+60+180))

# Copy the 300 sampled images to a folder to inspect them
file.copy(list3periods_shuffled,
          paste0("./Phase_1_check_sel_imgs/",
                 sapply(list3periods_shuffled, function(x) str_extract(x, "(SiteJE\\w+)")))
)

# Save the list of the 300 sampled images
write.csv(list3periods_shuffled,file = "./Phase_1_imgslist300.csv")

```
## Image patches labelling
For each image, a 200 pixels × 200 pixels image patch was randomly selected and plotted in RGB colours using the “plotRGB()” function of the “raster” package (Hijmans, 2022). Around 30 pixels per image were labelled by clicking on the image to retrieve the x and y coordinates using the “locator()” function of the “graphics” package and assigning to each pixel the class to which it belongs (see subsection 2.1). 

```r

library(raster)
library(rgdal)
rm(list=ls());
setwd("your/folder/path")

# load the image list and create seeds list that you will use to select a random 200X200 detail in the image
imgslist<-read.csv(file = "./Phase_1_imgslist300.csv",
                      row.names=1)[,1]
set.seed(1002); seedsx<-sample(1:2^15,300)
set.seed(1056); seedsy<-sample(1:2^15,300)
buffer<-100

#  PIXEL LABELLING PROCEDURE
#  1) define the classes
classes <- c("Soil", "Green_vegetation", "Kna_arv_flower", "Ran_acr_flower", "Leu_vul_flower", "Gra_flower")
#  2) select one class (e.g, Soil)
class<-classes[1]
#  3) run the line below where "RUN HERE" is written, an image will appear in the plot pane.
#  4) If there are some, click on some pixels of the class selected in point 2, (e.g., Soil),
     # then press esc on the keyboard.
     # A dataframe with coordinates, image name and label will be saved as a ".csv" file.
#  5) press Ctrl+Alt+T to run the section again. A new image will appear in the plot pane.
     # Repeat point 4) and 5) until all the images have been labelled
#  6) when you will have put some labels on all images, go to point 2) and select the second class.
     # Then go on with point 3),4),5),6) until all classes will be labelled.  
#  NB: LOCATOR WORKS PROPERLY ONLY WHEN ZOOM IN GLOBAL OPTIONS IS SET TO 100 %
   

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
                          labelled$type<-class;
                          labelled$im<-imgslist[i];
                          labelled$zoiCx<-zoiCx;
                          labelled$zoiCy<-zoiCy
                          write.csv(labelled,file=paste0("./Phase_1_lab_xy/","lab_xy_pic_",
                          sprintf("%03d", i),"class_",class,".csv"))}
                        }

```

## Close pixels removal
To prevent duplicated pixels after downscaling the images, any labelled pixels that were within a distance of 8 pixels from one another were eliminated from the dataset. The labelling phase resulted in a table where the class and pixel coordinates were stored.
```r
library(raster)
library(sp)
maindir<-"your/folder/path/"
setwd(maindir)
# combine all data frames into a single data frame
file_paths<-list.files(path="./Phase_1_lab_xy/",full.names=T,pattern=".csv")
combined_df <- do.call(rbind, lapply(file_paths, read.csv,row.names=1)
images_ids<-unique(combined_df$im)
minimum_distance<-8

# for each image load a first point, and then add only points more than 8 pixels apart
for(image in images_ids){
   df<-combined_df[combined_df$im==image,]
   ProxFiltered<- df[1,]   #at the first iteration only the first point is proximity filtered
   for (point_index in 2:nrow(df  )){
      pts<- as.matrix(  ProxFiltered[,c("x","y")])
      pt<-  as.numeric(df[point_index,c("x","y")])
      dists<- spDistsN1(pts, pt, longlat = FALSE)
      exceed<- any(dists<minimum_distance)
      if (exceed==FALSE){ProxFiltered<-rbind(  ProxFiltered, df[point_index,])}}
      if(!exists("ProxFiltered_tot")) {ProxFiltered_tot <- ProxFiltered} else {
      ProxFiltered_tot<-rbind(ProxFiltered_tot,ProxFiltered)};  
 }
# Save the labelled and proximity filtered dataset
write.csv(ProxFiltered_tot,"./Phase_1_labelled.csv")
```



# 2.1 Feature values computation
To increase the spectral separability of pixels between different classes, we computed RGB-based features: vegetation indices and texture metrics, described in detail in Table 1. We selected four vegetation indices well established in colour analysis literature (Lussem et al, 2018; Zhao, 2021). For pixels with a specific shade of purple colour, the calculation of the Visible Atmospherically Resistant Index (VARI) resulted in infinite values (for definition, see Table 1). Since only finite values can be used for classifier development, infinite VARI values were replaced with the highest finite value sampled (or lowest in case of negative infinite values), which occurred in less than 0.1% of the labelled pixels. The image textures were derived from co-occurrence matrices for each colour band, since we expected that the flower colours differed from the background (green vegetation or soil) surfaces (Guru et al., 2010), using the “glcm” package in R Studio (Zvoleff, 2020; Haralick et al, 1973). Homogeneity, contrast, dissimilarity, entropy, second moment, mean, and variance were computed in four directions (0°, 45°, 90° and 135°) and then averaged to one rotation-invariant texture as commonly used in texture analysis (e.g., Guru et al., 2010). For the computation of texture metrics, we needed to define the size of the window used for co-occurrence matrices. Moreover, downscaling the images to a lower resolution before feature extraction can give the best detection accuracy while also vastly increasing processing speed compared to higher resolution images (Mann et al., 2022).

We tested the influence of the size of the window used for co-occurrence matrices computation (WS) and downscaling factor (DF) on classification accuracy and speed. The tested values for WS were 3, 5, 7, 11, 19, 27, and 43 pixels (which is 3, 3+2, 3+4, 3+8, 3+16, 3+24, 3+40) and for DF, we used 1, 2, 4, and 8 pixel. Since flower dimensions were generally much smaller than one tenth of the image height, DF-WS combinations resulting in a window larger than one tenth of the image height were discarded (DF*WS<1040/10). The resulting number of DF-WS combinations was 23.

Here, we compute the 28 feature values for each of the downscaling factor - window size combinations. Prediction capability of random forest models developed based on the extracted features will be compared in phase 2.2
```r
library(dplyr)
library(randomForest)
library(crfsuite)
library(raster)
library(terra)
library(glcm)
setwd(maindir<-"your/folder/path")

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

setwd(maindir<-"your/folder/path")
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
load(paste("./ROISREFS/010/ROI/roi.data.Rdata",sep=""))


#------------------------------------------------------------
#  1 PROCESSING TIME OF RF MODELS INCLUDING ALL BANDS
#------------------------------------------------------------
  classifier_RF<-imgr<-imgm<-glcm.red<-glcm.green<-glcm.blue<-classified<-duration<-bestBands<-NULL
  classifier_RF<-readRDS("./Phase_3_RF_classifiers/RF_all.rds")
  time0 <- proc.time()[3]
    imgr<-brick("./IMGS/010/SiteJEA010_201404121334.jpg")
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


# 4.1 Flower covers extraction
Once the final RF classifier had been trained, the percentage of pixels in each class was computed for each image. For this, images were selected from the series (see section 2.3.1), for each image the selected features were computed, and percentages of each class within each image was calculated using the RF classifier developed in subsection 2.3.3. 

```r
library(randomForest)
library(glcm)
library(raster)
library(terra)
      setwd("your/folder/path")

# Define the feature extraction parameter, import the RF classifier, define plot IDs
    rf <- readRDS("./Phase_3_RF_classifiers/RF_selected.rds")
    ws=11
    df=4
    bestBands<-names(rf[["forest"]][["xlevels"]])
    bestBands
    plotIDs<-c( "001","002","003","005","007","008","009","010","011","013","016","017","018","021",
    "024","025","026","028","030","035","037","039","040","042","043","044","045","046","048","049",
    "051","053","054","056","057","058","059","060","062","064","065","067","070","071","073","074",
    "075","077","079","080","081","082","083","084","085","088","090","091","092","093","094","095",
    "097","099","100","102","103","105","108","109","110","111","113","114","115","116","119","120",
    "121","125","128","129","130","131","133","135","136","137","138")

# Extract the flower covers 
    # Expected processing time= 9 sec for each image, around 137 images per plot. 25 min/plot
    for (plot in plotIDs){
      df_sel<-read.csv("./Phase_1_selected_images.csv",row.names=1)
      imlist<-df_sel$imlistpath
      for (i in 1:length(imlist)){
          # Compute image features
            img<-brick(imlist[i])
            load(paste("./ROISREFS/",plot,"/ROI/roi.data.Rdata",sep=""))          
            imgr <- terra::aggregate(img, df)
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
    
          # Classify and compute the flower covers percentage on the not NA pixels (Areas not in the ROIs have NA values) 
            classified<-terra::predict(imgm,rf,na.rm=T)
            vals0<-summary(as.factor(classified@data@values))
            notNAcount<-sum(subset(vals0,!names(vals0)=="NA's"))       
            df_sel[i,rf[["classes"]][1]]<-if(length(subset(vals0,names(vals0)=="1"))>0){
                                                subset(vals0,names(vals0)=="1")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][2]]<-if(length(subset(vals0,names(vals0)=="2"))>0){
                                                subset(vals0,names(vals0)=="2")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][3]]<-if(length(subset(vals0,names(vals0)=="3"))>0){
                                                subset(vals0,names(vals0)=="3")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][4]]<-if(length(subset(vals0,names(vals0)=="4"))>0){
                                                subset(vals0,names(vals0)=="4")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][5]]<-if(length(subset(vals0,names(vals0)=="5"))>0){
                                                subset(vals0,names(vals0)=="5")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][6] ]<-if(length(subset(vals0,names(vals0)=="6"))>0){
                                                 subset(vals0,names(vals0)=="6")/notNAcount}else{NA}
      }
    # Save the flower covers of all selected images for each plot as a csv file.
    write.csv(df_sel,file = paste("./Phase_4_FCTS/dataflowers",plot,".csv",sep=""))
    }

```


# 4.2 FCTS smoothing
We identified and removed outliers from the derived flower cover time series using the “tsclean” function of the “forecast” R package which is based on Friedman's SuperSmoother for non-seasonal series (Hyndman & Khandakar, 2008). Values were aggregated at daily temporal resolution by by taking the arithmetic mean. A local polynomial regression function was fitted to smooth the time series using the “loess” function of the “stats” package (R Core Team, 2023).

```r
library(ggplot2)
library(zoo)
library(forecast)
library(dplyr)

setwd("your/folder/path")
list<-list.files(path="./Phase_4_FCTS/",pattern="dataflowers",full.names=T)
classes<-c("Gra_flower","Kna_arv_flower","Leu_vul_flower","Ran_acr_flower","Green_vegetation","Soil")
# Loop over each plot 
for (i in 1:length(list)){

  # Load the time-series
  plot<-substr(list[1],nchar(list[1])-6,nchar(list[1])-4)
  ts<-read.csv(list[i],stringsAsFactors = F)
  ts[is.na(ts)]<-0  #NA was assigned during classification to classes with zero pixels
  tszo<-read.zoo(ts[,c("date",classes)],tryFormats= c("%d/%m/%Y %H:%M","%Y-%m-%d %H:%M"))
  
   # Remove outliers in each class
     tot_df<-NULL
     for(c in 1:length(classes)){
            class<-classes[c]
            temp_df<-NULL
            temp_df<-fortify(tsclean(tszo[,classes[c]],replace.missing=F))
            temp_df$class<-class
            colnames(temp_df)[2]<-"FPe"
            if(is.null("tot_df")){tot_df<-temp_df}else{tot_df<-rbind(tot_df,temp_df)}
     }
     tot_df$date<-as.Date(tot_df$Index)
  
   # Aggregate to daily resolution
     tot_daily<-tot_df%>%group_by(date,class) %>%summarise_at(vars("FPe"), mean)
     
   # Add fitted values and standard error
     fitted_df<-NULL
     for(class in classes){
        temp_df<-NULL
        temp_df<-tot_daily[tot_daily$class==class,]
        if(nrow(temp_df)>0){
        temp_df<-cbind(temp_df,predict(loess(temp_df$FPe~as.numeric(as.Date(temp_df$date))),se=T))
              if(is.null("fitted_df")){fitted_df<-temp_df}else{fitted_df<-rbind(fitted_df,temp_df)}
        }
      }
     fitted_df$plot<-plot
     if(!exists("fall")){fall<-fitted_df}else{fall<-rbind(fitted_df,fall)}
}
   # Compute lower and upper confidence interval from fitted and standard error
      fall$lower<-fall$fit-qt(0.975,fall$df)*fall$se.fit
      fall$upper<-fall$fit+qt(0.975,fall$df)*fall$se.fit

write.csv(fall,file="./Phase_4_fitted_flower_cover_long_with_SE.csv")
```
# 4.3 Phenological metric identification
For each FCTS, onset, peak and end of flowering were extracted. The peak was identified as the day of maximum in the FCTS, where the value was higher than the values before and after it. The onset of flowering was identified on the basis of the rescaled cumulative sum of daily flower covers before the peak, whereas the end of flowering was identified on the basis of the rescaled cumulative sum of daily flower covers after the peak. This allowed the identification of flowering onset in FCTS when the end of the flowering was not observable (e.g., because of mowing) and the identification of the end of flowering in FCTS when the onset of flowering was not observable (e.g. image acquisition started later). The cumulative sums were scaled so that they ranged from 0 to 1. Onset was defined as the first day above 0.1 of the rescaled cumulative sum of daily flower covers before the peak, whereas the end of the season was identified as the first day above 0.9 of the rescaled cumulative sum of daily flower covers after the peak, as explained in Figure 3. The onset and the end of flowering were defined only in FCTS having a low flower cover (i.e., <1%) at the start and at the end of the observation period, respectively. All metrics were extracted only from time series where the peak was higher than 1%.

<figure>
<img src="https://drive.google.com/uc?id=1mPElue8oWmRnwIJ7BoIU2FAmIIoRke4z" width="600">
</figure>
            
_Flowering phenological metrics identification (Figure 4 in the manuscript)_
```r
library(dplyr)
setwd("your/folder/path")
classes<-c("Gra_flower","Kna_arv_flower","Leu_vul_flower","Ran_acr_flower","Green_vegetation","Soil")
classes_flowers<-classes[1:4]

# Upload all the fitted flower covers
fall<-read.csv(stringsAsFactors=T,file="./Phase_4_fitted_flower_cover_long_with_SE.csv",row.names=1)
fall$date<-as.Date(fall$date)
# keep only FCTS of flowers (remove "Green vegetation" and "Soil" covers)
fall<-fall[fall$class %in% classes_flowers,]

# Upload the plot composition, i.e. the list of species sown in each plot.
compo<-read.csv("./Phase_4_plot_composition.csv")
#create new column for each flower class, specifying if that flower species was sown in that plot
compo$Ran_acr_flower<-grepl("Ranunculus", compo$sp_list)
compo$Leu_vul_flower<-grepl("Leucanthemum", compo$sp_list)
compo$Kna_arv_flower<-grepl("Knautia", compo$sp_list)
compo$Gra_flower<-grepl("Festuca|Avenula|Poa|Anthoxanthum|Phleum|Dactylis|Holcus", compo$sp_list)

# Merge fitted FCTS and sown species list based on plot_ID
merged_df <- merge(fall, compo, by = "plot")

# Create a new column indicating if the flower species in the "class" column was sown in that plot
class_index <- match(merged_df$class, classes_flower)
merged_df$sown <- merged_df[, classes_flower][cbind(1:nrow(merged_df), class_index)]


limitonset<-.1
limitend<-.9

findpm<-function(df){
   onset<-peak<-end<-NA
   dates<-df$date
   vec<-df$fit
      dfbeforepeak=data.frame(date=dates[1:which.max(vec)],vec=vec[1:which.max(vec)])
      dfbeforepeak$rescaled<-rescale(cumsum(dfbeforepeak$vec))
      dfafterpeak<-data.frame(date=dates[(which.max(vec)):length(vec)],vec=vec[(which.max(vec)):length(vec)])
      dfafterpeak$rescaled<-rescale(cumsum(dfafterpeak$vec))
         onset<-dfbeforepeak$date[which(dfbeforepeak$rescaled>limitonset)[1]]
         peak<-dates[which.max(vec)]
         end<-dfafterpeak$date[which(dfafterpeak$rescaled>limitend)[1]]
         pm<-c(onset,peak,end)
           if(max(vec)<0.01){pm<-NA}
           if(vec[1]>(.01)){pm[1]<-NA}
           if(vec[length(vec)]>(.01)){pm[3]<-NA}
           if(peak==max(dates)|peak==min(dates)){pm<-NA}
   return(pm)
}


# Identify flowering phenological metric for each FCTS 
phenometrics<-merged_df%>%
   group_by(plot,class,sown) %>%
   dplyr::summarize(onset = findpm(cur_data())[1],
                    peak = findpm(cur_data())[2],
                    end = findpm(cur_data())[3])

# Identify flowering phenological metric for FCTS of the sown species
phenometrics_sown<-phenometrics%>%
   filter(sown==T)
summary(is.na(phenometrics_sown))
```


# IMAGE CLASSIFICATION
Here, we will classify and save the images of one plot. 
![Example of RGB and classified image from Figure 5 in the manuscript](https://drive.google.com/uc?id=1NVcvDAzGoqVIJ4gtlL2xgXSbHAvd22ZY)  
*Example of RGB and classified image from Figure 5 in the manuscript*

```r
library(randomForest)
library(raster)
library(glcm)
library(terra)
library(stringr)
setwd("your/folder/path/")

# Upload the RF classifier and define the plot you want to classify
# and create a folder where the classified images will be stored
rf <- readRDS("./Phase_3_RF_classifiers/RF_selected.rds")
plot<-"056"
dir<-paste0("./Phase_4_Classified_images/",plot)
dir.create(dir)

# Read the CSV containing the selected images IDs based on brightness and contrast
df_sel_im_all<-read.csv("./Phase_1_selected_images.csv")
pathimlist<-df_sel_im_all[substr(df_sel_im_all$imname,8,10)==plot,"imlistpath"]

# Define the ws and df, upload the images, compute the features,
# classify the image and save a false color plot of the classified image
ws=11
df=4

for (imgn in pathimlist){
   img<-brick(imgn)
   imgr <- terra::aggregate(img, df)
   #NB: Here compute only the "best Bands selected in SFFS" in our case:
   #    rgbvi,gli,vari,ngrdi, R_second_moment, B_contrast, B_entropy and B_second_moment
      imgr[[4]]<-((imgr[[2]]*imgr[[2]])-(imgr[[1]]*imgr[[3]]))/((imgr[[2]]*imgr[[2]])+(imgr[[1]]*imgr[[3]]))#RGBVI
      imgr[[5]]<-((2*imgr[[2]])-imgr[[1]]-imgr[[3]])/((2*imgr[[2]])+imgr[[1]]+imgr[[3]])#GLI
      imgr[[6]]<-(imgr[[2]]-imgr[[1]])/(imgr[[2]]+imgr[[1]]-imgr[[3]])#VARI
      imgr[[6]][imgr[[6]]==Inf]<-max(getValues(imgr[[6]])[is.finite(getValues(imgr[[6]]))])
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
      rasterVis::levelplot(classified,  margin=F,colorkey=F,scales=list(draw=F),
                at = seq(1,6), col.regions=c("grey20","seagreen","purple","cyan","gold","brown"))
   dev.off()
   #----------------Copy the original images in the same subfolder to compare classified and RGB
   file.copy(from=imgn,to=basename(imgn))
}
```


