# 1.2 Images selection based on brightness and contrast
Images with uniform light conditions were retrieved by selecting brightness and contrast between the 10th and the 40th percentile within a 3-day window. The selection of the best images within this 3-day window avoided including images taken on days with sub-optimal observations (e.g., all foggy/high contrast images).

```r
#####----------- CREATE LIST OF BLOCKS AND PLOTS ----------###
plotlist<-c(
   c("001", "002", "003", "005", "007", "008", "009", "010", "011", "013", "016", "017", "018", "020", "021", "024", "025", "026", "027", "028", "030", "033", "035", "037", "039", "040", "042", "043", "044" ,"045", "046"),
   c("048", "049", "051", "053", "054", "056", "057", "058", "059", "060", "062", "064", "065", "067", "070", "071", "073", "074", "075", "077", "079", "080", "081", "082", "083", "084", "085", "088", "090", "091", "092"),
   c("093", "094", "095", "097", "099", "100", "102", "103", "105", "108", "109", "110", "111", "113", "114", "115", "116", "119", "120", "121", "125", "128", "129", "130", "131", "133", "135", "136", "137", "138")
     )
blocklist<- c(rep("A",31), rep("B",31) ,rep("C",30))
blockplotlist<-paste(blocklist,plotlist,sep="")
#now remove plots where the image acquisition failed
blockplotlist<-blockplotlist[!blockplotlist %in% c("A020","A027","A033")]
count<-data.frame(blockplotlist,NA,NA)
colnames(count)<-c("bp","totN","selectedN")


start=113           #Start of the season of interest (doy)
end=150             #End of the season of interest (doy)
perc<-c(0.1,0.4)    #Images with uniform light conditions were retrieved by selecting brightness and contrast between the 10th and the 40th percentile in a 3-days window

#######------LOOP OVER THE SELECTED PLOTS -------####
for (b in 1:length(blockplotlist)){
  bp<-blockplotlist[b]
  block<- substr(bp,1,1)
  plot<- substr(bp,2,4)
   dataraw<- read.csv(paste("your/folder/path/Phase_1_2014_raw_indices/VI raw2014 ",block,"-",plot,".csv",sep=""),
                       encoding="UTF-8",row.names=1 )
   doyslist<-dataraw$doy[!duplicated(dataraw$doy)]
   doyslist<-doyslist[doyslist>start&doyslist<end]
   # Prepare empty dataframe
      dataselected0<-dataselectedsd0<-dataraw[0,]

#######------FIND 0.1 to 0.4 percentile bri.av in 3 days window-------####
for (i in 1:length(doyslist)){
    tab3days<-dataraw[dataraw$doy>doyslist[i]-2&dataraw$doy<doyslist[i]+2,]
    tab3days$percentiles<-ecdf(tab3days$bri.av)(tab3days$bri.av)
    dataselected0<-rbind(dataselected0,
                         tab3days[tab3days$percentiles>perc[1]&tab3days$percentiles<perc[2],])

}
#eliminate duplicated
datafs<-dataselected0[!duplicated(dataselected0$date), ]

#ADD a column identifying images matching the criteria of selection for bri.av
match<-is.element(dataraw$date, datafs$date)
sum(match)
dataraw$selBRIAV<-0
dataraw$selBRIAV[match]<-1

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
