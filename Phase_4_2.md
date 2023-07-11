# 4.2 FCTS smoothing and display
## FCTS smoothing
We identified and removed outliers from the derived flower cover time series using the “tsclean” function of the “forecast” R package which is based on Friedman's SuperSmoother for non-seasonal series (Hyndman & Khandakar, 2008). Values were aggregated at daily temporal resolution by averaging. A local polynomial regression function was fitted to smooth the time series using the “loess” function of the “stats” package (R Core Team, 2023).

```r

library(ggplot2)
library(zoo)
library(forecast)
library(reshape2)
library(dplyr)
library(ggtext)
library(ggpubr)
Sys.setlocale("LC_TIME", "English")

#-----------------------------------------------------------------------------
#---------------------------1) FROM RAW TIME-SERIES TO FITTED TIME-SERIES
#-----------------------------------------------------------------------------
rm(list=ls())
setwd("your/folder/path/Phase_4_FCTS/")
list<-list.files(pattern="dataflowers",full=F)

classes<-c("Gra_flower","Kna_arv_flower","Leu_vul_flower","Ran_acr_flower","Green_vegetation","Soil")
for (i in 1:length(list)){

 #---------------------------1.1) LOAD THE TIME_SERIES and REMOVE OUTLIERS
    # "Gra_flower","Kna_arv_flower"+"Leu_vul_flower"+"Ran_acr_flower"
  plot<-substr(list[i],12,14)
  plot
  setwd("your/folder/path/Phase_4_FCTS/")
  ts<-read.csv(list[i],stringsAsFactors = F)
  ts[is.na(ts)]<-0  #NA was assigned during classification to classes with zero pixels in the image, so we can set it to zero
  tszo<-read.zoo(ts[,c("date",classes)],tryFormats= c("%d/%m/%Y %H:%M","%Y-%m-%d %H:%M"))
  
   # remove outliers
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
     
     
     #-----------------------1.2) add fitted values and standard error
     fitted_df<-NULL
      for(c in 1:length(classes)){
        class<-classes[c]
        temp_df<-NULL
        temp_df<-tot_daily[tot_daily$class==class,]
        if(nrow(temp_df)>0){
        temp_df<-cbind(temp_df,predict(loess(temp_df$FPe~as.numeric(as.Date(temp_df$date))),se=T))
        if(is.null("fitted_df")){fitted_df<-temp_df}else{fitted_df<-rbind(fitted_df,temp_df)}
     }}
   
    fitted_df$plot<-plot
  if(!exists("fall")){fall<-fitted_df}else{fall<-rbind(fitted_df,fall)}
    }
   #Compute lower and upper confidence interval from fitted and standard error
      fall$lower<-fall$fit-qt(0.975,fall$df)*fall$se.fit
      fall$upper<-fall$fit+qt(0.975,fall$df)*fall$se.fit

write.csv(fall,file="your/folder/path/Phase_4_fitted_flower_cover_long_with_SE.csv")

```
## FCTS smoothing
The following code creates a figure of the FCTS for each of the 89 plots

```
library(ggplot2)
library(dplyr)
setwd("your/folder/path/Phase_4_FCTS_plots/")

rm(list=ls())
classes<-c("Gra_flower","Kna_arv_flower","Leu_vul_flower","Ran_acr_flower","Green_vegetation","Soil")
classes_display<-classes[1:4]
colors<-c("grey20","purple","cyan","gold","forestgreen","brown")
fall0<-read.csv(stringsAsFactors=T,file="your/folder/path/Phase_4_fitted_flower_cover_long_with_SE.csv",row.names=1)
fall0$date<-as.Date(fall0$date)
fall0<-fall0[fall0$class %in% classes_display,]


for (i in unique(fall0$plot)){
   fall<-fall0[fall0$plot==i,]
   resi<-res[res$plot==i,]
p<-ggplot(data=fall,aes(date,fit,fill=factor(class,level=classes),colour=factor(class,level=classes)))+
  ggtitle(paste0("Plot ",fall$plot[1]))+
  ylab("Flower cover")+
  geom_line()+
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,col="transparent")+
  coord_cartesian(ylim=c(0, if(max(fall$FPe)>0.3){max(fall$FPe)}else{0.3})) +
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=colors,#labels=classes
                       labels=c("Graminoids","*K. arvensis*", "*L. vulgare*","*R. acris*"))+
  scale_color_manual(values=colors,#labels=classes
                     labels=c("Graminoids","*K. arvensis*", "*L. vulgare*","*R. acris*"))+
  theme_bw()+
  guides(fill=guide_legend(nrow=2,byrow=TRUE,title="Legend"),
         color=guide_legend(nrow=2,byrow=TRUE,title="Legend"))+
  theme(legend.position="bottom",
        legend.background = element_rect(fill="transparent"),
        panel.grid = element_blank(),
        legend.text = element_markdown(size=10))

p
ggsave(plot=p, filename=paste0("FCTS_plot",i,".png"),
       units="cm", width=12, height=12,dpi=320)
}

