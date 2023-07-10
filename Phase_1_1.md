# 1.1 Brightness and contrast extraction

Light conditions heavily affect pixel colours: images with high brightness were usually foggy, and images with high contrast were usually acquired in direct sunlight conditions. We calculated brightness and contrast for all images using the “extractVIs” function of the R package “Phenopix” (Filippa et al., 2016) and tested which brightness and contrast combinations allow the selection of images acquired in homogeneous light conditions.
Directory structure in the "your/path/folder/ROISREFS_2014" is compatible with the one required in phenopix, a very commonly used package for phenocam images processing. The extraction of each time-series required around 20 minutes.

In this script, our aim is to extract brightness average and contrast for each image.

```
library(phenopix)
#####----------------------------------------------------------------------###
#####----------- 1) create list of paths to images in subfolders ----------###
#####----------------------------------------------------------------------###

  plots.A <- c("001", "002", "003", "005", "007", "008", "009", "010", "011", "013", "016", "017", "018", "020", "021", "024", "025", "026", "027", "028", "030", "033", "035", "037", "039", "040", "042", "043", "044" ,"045", "046")
  plots.B <- c("048", "049", "051", "053", "054", "056", "057", "058", "059", "060", "062", "064", "065", "067", "070", "071", "073", "074", "075", "077", "079", "080", "081", "082", "083", "084", "085", "088", "090", "091", "092")
  plots.C <- c("093", "094", "095", "097", "099", "100", "102", "103", "105", "108", "109", "110", "111", "113", "114", "115", "116", "119", "120", "121", "125", "128", "129", "130", "131", "133", "135", "136", "137", "138")

  plots.A.2 <- paste("A", plots.A, sep="")
  plots.B.2 <- paste("B", plots.B, sep="")
  plots.C.2 <- paste("C", plots.C, sep="")

  all.plots <- c(plots.A.2, plots.B.2, plots.C.2)
  all.plots.path <- paste(substr(all.plots, 1, 1), "/", substr(all.plots, 2, 4), sep="")

#####----------------------------------------------------------------------###
#####----------- 2) Extract Vegetation Indices-----------------------------###
#####---------------around 20 min per plot---------------------------------###
  
    time0.all <- proc.time()[3]

    for (plot in all.plots.path) {
                time0 <- proc.time()[3]
                path.JPG <- paste("your/folder/path/IMGS_2014/", substr(plot, 3,5), "/", sep="")
                path.ROI <- paste("your/folder/path/ROISREFS_2014/", plot, "/ROI/", sep="")
                path.VI <- paste("your/folder/path/ROISREFS_2014/", plot, "/VI/", sep="")
                setwd(path.JPG); print(getwd())
               extracted<-extractVIs(path.JPG, path.ROI, vi.path=path.VI,plot=F, spatial=F, date.code="yyyymmddHHMM")
               duration <- proc.time()[3]-time0;                print(duration)

         }

time1.all <- proc.time()[3]
duration.all <- time1.all-time0.all 



#####----------------------------------------------------------------------###
#####----------- from Rdata format to csv format---------------------------###
#####----------------------------------------------------------------------###

       for (plot in all.plots.path) {
      load(paste0("your/folder/path/ROISREFS_2014/", plot, "/VI/VI.data.Rdata"))
      path.csv <- "your/folder/path/Phase_1_2014_raw_indices/"
      write.csv(VI.data$roi1,file=paste0(path.csv,"VI raw2014 ", substr(plot, 1,1), "-",substr(plot, 3,5), ".csv"))
    }
```
