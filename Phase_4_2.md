# 4.2 FCTS smoothing and display
We identified and removed outliers from the derived flower cover time series using the “tsclean” function of the “forecast” R package which is based on Friedman's SuperSmoother for non-seasonal series (Hyndman & Khandakar, 2008). Values were aggregated at daily temporal resolution by averaging. A local polynomial regression function was fitted to smooth the time series using the “loess” function of the “stats” package (R Core Team, 2023).
```r
```
