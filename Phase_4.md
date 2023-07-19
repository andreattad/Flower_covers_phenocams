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
    plotIDs<-c("001", "002", "003", "005", "007", "008", "009", "010", "011", "013", "016", "017", "018",  "021", "024", "025", "026", "028", "030", "035", "037", "039", "040", "042", "043", "044" ,"045", "046","048", "049", "051", "053", "054", "056", "057", "058", "059", "060", "062", "064", "065", "067", "070", "071", "073", "074", "075", "077", "079", "080", "081", "082", "083", "084", "085", "088", "090", "091", "092","093", "094", "095", "097", "099", "100", "102", "103", "105", "108", "109", "110", "111", "113", "114", "115", "116", "119", "120", "121", "125", "128", "129", "130", "131", "133", "135", "136", "137", "138")

# Extract the flower covers 
    # Expected processing time= 9 sec for each image, around 137 images per plot. Around 25 minutes for each plot
    for (plot in plotIDs){
      df_sel$Gra_flower<-df_sel$Green_vegetation<-df_sel$Kna_arv_flower<-df_sel$Leu_vul_flower<-df_sel$Ran_acr_flower<-df_sel$Soil<-NA
      time0 <- proc.time()[3]; 
      df_sel<-read.csv("./Phase_1_selected_images.csv",row.names=1)
      imlist<-df_sel$imlistpath
      for (i in 1:length(imlist)){
          # Compute image features
            img<-brick(imlist[i])
            load(paste("./ROISREFS_2014/",plot,"/ROI/roi.data.Rdata",sep=""))          
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
    
          # Classify and compute the flower covers percentage on the not NA pixels (Areas not in the ROIs have NA values) 
            classified<-terra::predict(imgm,rf,na.rm=T)
            vals0<-summary(as.factor(classified@data@values))
            notNAcount<-sum(subset(vals0,!names(vals0)=="NA's"))       
            df_sel[i,rf[["classes"]][1]]<-if(length(subset(vals0,names(vals0)=="1"))>0){subset(vals0,names(vals0)=="1")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][2]]<-if(length(subset(vals0,names(vals0)=="2"))>0){subset(vals0,names(vals0)=="2")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][3]]<-if(length(subset(vals0,names(vals0)=="3"))>0){subset(vals0,names(vals0)=="3")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][4]]<-if(length(subset(vals0,names(vals0)=="4"))>0){subset(vals0,names(vals0)=="4")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][5]]<-if(length(subset(vals0,names(vals0)=="5"))>0){subset(vals0,names(vals0)=="5")/notNAcount}else{NA}
            df_sel[i,rf[["classes"]][6] ]<-if(length(subset(vals0,names(vals0)=="6"))>0){subset(vals0,names(vals0)=="6")/notNAcount}else{NA}
            time1 <- proc.time()[3]
            duration <- time1-time0
            print(duration)
      }
    # Save the flower covers of all selected images for each plot as a csv file.
    write.csv(df_sel,file = paste("./Phase_4_FCTS/dataflowers",plot,".csv",sep=""))
    }

```


# 4.2 FCTS smoothing
We identified and removed outliers from the derived flower cover time series using the “tsclean” function of the “forecast” R package which is based on Friedman's SuperSmoother for non-seasonal series (Hyndman & Khandakar, 2008). Values were aggregated at daily temporal resolution by averaging. A local polynomial regression function was fitted to smooth the time series using the “loess” function of the “stats” package (R Core Team, 2023).

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

  # LOAD THE TIME_SERIES and REMOVE OUTLIERS
  plot<-substr(list[1],nchar(list[1])-6,nchar(list[1])-4)
  ts<-read.csv(list[i],stringsAsFactors = F)
  ts[is.na(ts)]<-0  #NA was assigned during classification to classes with zero pixels in the image, so we can set it to zero
  tszo<-read.zoo(ts[,c("date",classes)],tryFormats= c("%d/%m/%Y %H:%M","%Y-%m-%d %H:%M"))
  
   # remove outliers in each class
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
  
   # aggregate to daily resolution
     tot_daily<-tot_df%>%group_by(date,class) %>%summarise_at(vars("FPe"), mean)
     
   # add fitted values and standard error
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
   #Compute lower and upper confidence interval from fitted and standard error
      fall$lower<-fall$fit-qt(0.975,fall$df)*fall$se.fit
      fall$upper<-fall$fit+qt(0.975,fall$df)*fall$se.fit

write.csv(fall,file="./Phase_4_fitted_flower_cover_long_with_SE.csv")
```
# 4.3 Phenological metric extraction
For each FCTS, onset, peak and end of flowering were extracted. The peak was identified as the day of maximum in the FCTS. Onset of flowering was identified as the first day above the 10th percentile of the cumulative flower cover until peak, whereas the end of the season was identified as the first day above the 90th percentile of the cumulative flower cover after the peak, as explained in Figure 3. In our case study it was not possible to extract all the metrics from all FCTS, since in a few cases flowering started before the start of observations (i.e. spring weeding), and in many cases the peak and the end of flowering were not observed because mowing interrupted grassland development. Therefore, we limited the flowering metrics extraction: the onset was extracted in FCTS having flower cover at the first observation day lower than one third of the maximum flower cover, end of the season in FCTS having flower cover at the last observation day lower than one third of the maximum flower cover, all metrics were extracted only in FCTS with peak after first observation day and before last observation day and flower cover at peak higher than 1%.

<figure>
  <img src="https://drive.google.com/uc?id=1mPElue8oWmRnwIJ7BoIU2FAmIIoRke4z" width="600">
  <figcaption> __Flowering phenological metrics identification__ </figcaption>
</figure>

```r
library(dplyr)
setwd("your/folder/path")
classes<-c("Gra_flower","Kna_arv_flower","Leu_vul_flower","Ran_acr_flower","Green_vegetation","Soil")
classes_flowers<-classes[1:4]

# Upload all the fitted flower covers
fall<-read.csv(stringsAsFactors=T,file="./Phase_4_fitted_flower_cover_long_with_SE.csv",row.names=1)
fall$date<-as.Date(fall0$date)
# keep only FCTS of flowers (remove "Green vegetation" and "Soil" covers)
fall<-fall[fall$class %in% classes_flowers,]

# Upload the plot composition, i.e. the list of species sown in each plot.
compo<-read.csv(your/folder/path/Phase_4_plot_composition.csv")
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
