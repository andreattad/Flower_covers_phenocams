# 4.3 Phenological metric extraction
For each FCTS, onset, peak and end of flowering were extracted. The peak was identified as the day of maximum in the FCTS. Onset of flowering was identified as the first day above the 10th percentile of the cumulative flower cover until peak, whereas the end of the season was identified as the first day above the 90th percentile of the cumulative flower cover after the peak, as explained in Figure 3. In our case study it was not possible to extract all the metrics from all FCTS, since in a few cases flowering started before the start of observations (i.e. spring weeding), and in many cases the peak and the end of flowering were not observed because mowing interrupted grassland development. Therefore, we limited the flowering metrics extraction: the onset was extracted in FCTS having flower cover at the first observation day lower than one third of the maximum flower cover, end of the season in FCTS having flower cover at the last observation day lower than one third of the maximum flower cover, all metrics were extracted only in FCTS with peak after first observation day and before last observation day and flower cover at peak higher than 1%.

![Flowering phenological metrics identification (Figure 3 in the manuscript)](pm_extraction.png)

```r
```
