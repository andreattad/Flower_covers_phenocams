# **Tutorial for the manuscript "Extracting flowering phenology from grassland species mixtures using time-lapse cameras"**
by Davide Andreatta, Christoph Bachofen, Michele Dalponte, Valentin H. Klaus, Nina Buchmann.   
Published in Remote Sensing of Environment    
DOI: https://doi.org/10.1016/j.rse.2023.113835   
Contact information: davide.andreatta@phd.unipd.it.    

<figure>
<img src="https://drive.google.com/uc?id=1ByJcwP-dgLWCyJVBRwHEPDjsj_EQEqUw" width="800">
</figure>

_Example of RBG images, classified images and flower cover time-series_




In this tutorial single species flowering phenology time-series and phenological metrics from time-lapse camera images of grasslands captured at the Jena trait-based experiment (Germany) were extracted. 

To ensure that the results of the manuscript can be easily replicated, we decided to base our tutorial on the provided case study and provide code that can be immediately applied to the supplementary materials, including images, labels, and regions of interest, which will be made available on the ETH repository (https://doi.org/10.3929/ethz-b-000634004). Although the presented workflow is tailored to the specific case study, it can be adapted for use with different case studies with minor modifications. To apply the workflow to new case studies, adjustments should be made to plot IDs, to the period of interest, and to flower class names.

The code is developed in R version 4.3.0 (2023-04-21 ucrt).

<figure>
<img src="https://drive.google.com/uc?id=1z9yc0Tht4S425txHhc42OpQdz4N7g6I6" width="800">
</figure>

_Structure of the proposed workflow (Figure 2 in the manuscript). Abbreviations are as follows: RF = Random Forest; SFFS = Sequential Floating Forward Selection; FCTS = Flower Cover Time Series._

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
The materials can be downloaded at the link: https://doi.org/10.3929/ethz-b-000634004     
Materials include: images, region of interests and labels.
1) Download in your working directory the 4 zip files available in the ETH REPOSITORY (https://doi.org/10.3929/ethz-b-000634004)     
2) Create a new folder and name it "path/to/your/working/directory/IMGS"
3) Unzip the zip folders "IMGS_PLOTS_001_049.zip", "IMGS_PLOTS_051_093.zip", "IMGS_PLOTS_094_138.zip"
4) Put the content of the unzipped folders in "path/to/your/working/directory/IMGS"
5) unzip "ROIS_LABELS_COMPOSITION.zip"

## 2. Set up the working environment. 

```
setwd("path/to/your/working/directory/")
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
