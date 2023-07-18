# 1.1 Brightness and contrast extraction

Light conditions heavily affect pixel colours: images with high brightness were usually foggy, and images with high contrast were usually acquired in direct sunlight conditions. We calculated brightness and contrast for all images using the “extractVIs” function of the R package “Phenopix” (Filippa et al., 2016) and tested which brightness and contrast combinations allow the selection of images acquired in homogeneous light conditions.
Directory structure in the "your/path/folder/ROISREFS_2014" is compatible with the one required in phenopix, a very commonly used package for phenocam images processing. The extraction of each time-series required around 20 minutes.

In this script, our aim is to extract brightness average and contrast for each image.

```r
library(phenopix)

plotIDs<-c("001", "002", "003", "005", "007", "008", "009", "010", "011", "013", "016", "017", "018",  "021", "024", "025", "026", "028", "030", "035", "037", "039", "040", "042", "043", "044" ,"045", "046","048", "049", "051", "053", "054", "056", "057", "058", "059", "060", "062", "064", "065", "067", "070", "071", "073", "074", "075", "077", "079", "080", "081", "082", "083", "084", "085", "088", "090", "091", "092","093", "094", "095", "097", "099", "100", "102", "103", "105", "108", "109", "110", "111", "113", "114", "115", "116", "119", "120", "121", "125", "128", "129", "130", "131", "133", "135", "136", "137", "138")

# Extract vegetation indices
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

# Convert from Rdata format to csv format and save all in the same folder
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

#######------FIND 0.1 to 0.4 percentile bri.av in 3 days window-------####
count<-data.frame(plotIDs=plotIDs,totN=NA,selectedN=NA)
for (plot in 1:length(plotIDs)){
   dataraw<- read.csv(paste0("your/folder/path/Phase_1_2014_raw_indices/"VI_raw2014_", plot,".csv"),encoding="UTF-8",row.names=1 )
   doys<-unique(dataraw$doy)
   doys<-doys[doys>start&doys<end]
   for (i in 1:length(doys)){
        tab3days<-dataraw[dataraw$doy>doys[i]-2&dataraw$doy<doys[i]+2,]
        tab3days$percentiles<-ecdf(tab3days$bri.av)(tab3days$bri.av)
        dates_selectedi<-tab3days[tab3days$percentiles>perc[1]&tab3days$percentiles<perc[2],])
        if(i==1){dates_selectedtot<-dates_selectedi}else{
                 dates_selectedtot<-rbind(dates_selectedtot,dates_selectedi)
        }
# eliminate duplicated
dates_selectedtot<-dates_selectedtot[!duplicated(dates_selectedtot$date), ]

# add column selBRI in which images matching the criteria of selection for brightness are TRUE 
match<-is.element(dataraw$date, dates_selectedtot$date)
dataraw$selBRI<-FALSE
dataraw$selBRI[match]<-TRUE

#######------FIND 0.1 to 0.4 percentile bri.sd in 3 days window-------####
for (i in 1:length(doyslist)){
   tab3dayssd<-dataraw[dataraw$doy>doyslist[i]-2&dataraw$doy<doyslist[i]+2,]
   tab3dayssd$percentiles<-ecdf(tab3dayssd$bri.sd)(tab3dayssd$bri.sd)
   dataselectedsd0<-rbind(dataselectedsd0,
                        tab3dayssd[tab3dayssd$percentiles>perc[1]&tab3dayssd$percentiles<perc[2],])
}
#eliminate duplicated
datafssd<-dataselectedsd0[!duplicated(dataselectedsd0$date), ]

#ADD a column identifying images matching the criteria of selection for bri.sd
matchsd<-is.element(dataraw$date, datafssd$date)
sum(match)
dataraw$selBRISD<-0
dataraw$selBRISD[matchsd]<-1


#######------ADD a COLUMN to identify images respecting both the conditions-------####

dataraw$selSUM<-dataraw$selBRISD+dataraw$selBRIAV
dataraw$selSUM[dataraw$selSUM==1]<-0
dataraw<-dataraw[dataraw$doy>start&dataraw$doy<end,]
dataraw$imname<-paste0("SiteJE",block,plot,"_2014",substr(dataraw$date,6,7),#yearmonth
                       substr(dataraw$date,9,10),substr(dataraw$date,12,13),#dayhour
                       substr(dataraw$date,15,16),".jpg")#minute
count[b,2]<-nrow(dataraw)
count[b,3]<-nrow(dataraw[dataraw$selSUM==2,])
 setwd("your/folder/path/Phase_1_2014_filtered_indices_BRIAV_BRISD/")
 write.csv(dataraw,file = paste("rawdatafilt",plot,".csv",sep=""))
}

#######------write a summary csv-------####
write.csv(count,file = "your/folder/path/Phase_1_Available_selected_images_per_plot.csv")

```

