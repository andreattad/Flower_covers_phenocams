# 1.1 Brightness and contrast extraction

Light conditions heavily affect pixel colours: images with high brightness were usually foggy, and images with high contrast were usually acquired in direct sunlight conditions. We calculated brightness and contrast for all images using the “extractVIs” function of the R package “Phenopix” (Filippa et al., 2016) and tested which brightness and contrast combinations allow the selection of images acquired in homogeneous light conditions.
Directory structure in the "your/path/folder/ROISREFS_2014" is compatible with the one required in phenopix, a very commonly used package for phenocam images processing. The extraction of each time-series required around 20 minutes.

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
  path.JPG <- paste("./IMGS_2014/", plot"/", sep="")
  path.ROI <- paste("./ROISREFS_2014/", plot, "/ROI/", sep="")
  path.VI <- paste("./ROISREFS_2014/", plot, "/VI/", sep="")
  extracted<-extractVIs(path.JPG, path.ROI, vi.path=path.VI,plot=F, spatial=F, date.code="yyyymmddHHMM")
  duration <- proc.time()[3]-time0;
  print(duration)
}

# Convert from Rdata format to csv format and save all them in the same folder
for (plot in plotIDs) {
  load(paste0("./ROISREFS_2014/", plot, "/VI/VI.data.Rdata"))
  path.csv <- "./Phase_1_2014_raw_indices/"
  write.csv(VI.data$roi1,file=paste0(path.csv,"VI_raw2014_", plot,".csv"))
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
   dataraw<- read.csv(paste0("./Phase_1_2014_raw_indices/VI_raw2014_", plot,".csv"),
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
               write.csv(dataraw,file = paste("./Phase_1_2014_filtered_indices_BRI_CON/rawdatafilt",plot,".csv",sep=""))
}

```

# 1.3 Image labelling
Image labelling phase includes: i) preparation of the list of images to be labelled, ii) label the images, iii) remove labelled pixels which are too close from one another. The labelled dataset is provided in the ETH repository (see README file for details) for replicability.

## Image list preparation
To develop a labelled dataset, 300 images were randomly selected (60 images in the period between Apr, 24 and May, 5; 60 images between May, 6 and May, 18; 120 images between May, 19 and May, 29). 
```r
library(stringr)

setwd("your/folder/path")

# Create a list of all the filtered images
csv_name_list<-list.files(path="./Phase_1_2014_filtered_indices_BRI_CON/",pattern="rawdatafilt")
combined.df <- do.call(rbind , lapply(csv_name_list, read.csv, row.names = 1))
combined.df_sel<-combined.df[combined.df$selBRICON==T,]
combined.df_sel$imlistpath<-paste0(maindir,"IMGS_2014/",
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
## Image patches  labelling
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
