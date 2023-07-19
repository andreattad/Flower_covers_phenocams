# **Tutorial for the manuscript "Extracting single species flowering phenology from grassland species mixtures using time-lapse cameras"**
by D.Andreatta, C. Bachofen, M. Dalponte, V. Klaus, N.Buchmann.   
Submitted to Remote Sensing of Environment. 
Contact information: davide.andreatta@phd.unipd.it.

In this tutorial single species flowering phenology time-series and phenological metrics from time-lapse camera images of grasslands captured at the Jena trait-based experiment (Germany) were extracted. 
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
     -    3.1 feature selection and models comparison
     -    3.2 processing time calculation
4. [FCTS extraction, smoothing and calculation of phenological metrics.](Phase_4.md)
     -    4.1 FCTS extraction
     -    4.2 FCTS smoothing and display
     -    4.3 Phenological metric extraction
  
Moreover, [here](image_classification.md) we provide the code to classify all the images in a plot
![Example of RGB and classified image from Figure 5 in the manuscript](https://drive.google.com/uc?id=1NVcvDAzGoqVIJ4gtlL2xgXSbHAvd22ZY)  
*Example of RGB and classified image*

# To get ready for the tutorial:

## 1. Download the materials:
the materials can be downloaded at the link: LINK TO ETH REPOSITORY HERE
materials include: Images, Region of interests and labels.

## 2. Set up the working environment. 

```
setwd("your/folder/path")
dir.create("Phase_1_2014_raw_indices")
dir.create("Phase_1_2014_filtered_indices_BRI_CON")
dir.create("Phase_1_lab_xy")
dir.create("Phase_2_df_ws_ext_feat")
dir.create("Phase_3_RF_classifiers")
dir.create("Phase_4_FCTS")
dir.create("Phase_4_Classified_images_")
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
