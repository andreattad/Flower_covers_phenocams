
# **TUTORIAL for the manuscript "Extracting single species flowering phenology from grassland species mixtures using time-lapse cameras"**
 by D.Andreatta, V. Klaus, C. Bachofen, M. Dalponte, N.Buchmann.  
 submitted to Remote Sensing of Environment

In this tutorial we will derive single species flowering phenology time-series and phenological metrics from time-lapse camera images of grasslands. 
The code is developed in R version 4.3.0 (2023-04-21 ucrt)


to get ready for the tutorial:




## 1. Download the materials:
the materials can be downloaded at the link: LINK TO ETH REPOSITORY HERE
materials include: Images, Region of interests and labels. Download them in "your/folder/path"

## 2. Install the required R packages

```

packages <- c("bigstatsr", "crfsuite", "data.table", "dplyr", "forecast", "ggplot2", 
              "ggpubr", "ggrepel", "ggtext", "glcm", "grid", "gridExtra", "jpeg",
              "lattice", "phenopix", "randomForest", "raster", "readr", "reshape2",
              "rgdal", "sp", "stringr", "terra", "tidyr", "varSel", "zoo")

# Load or install missing packages
lapply(packages, function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
})

```

## 3. Set up the working environment. 

```
setwd("your/folder/path")
dir.create("Phase_1_2014_raw_indices")
dir.create("Phase_1_2014_filtered_indices_BRIAV_BRISD")
dir.create("Phase_1_lab_xy")
dir.create("Phase_2_df_ws_ext_feat")
dir.create("Phase_3_RF_classifiers")
dir.create("Phase_4_2014_filtered_indices_GCC")
dir.create("Phase_4_Classified_images_for_Fig_5")
dir.create("Phase_4_FCTS")
dir.create("Phase_4_FCTS_plots")

```
