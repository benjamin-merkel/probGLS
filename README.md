# probGLS

The package probabilistic algorithm for geolocation data (probGLS) provides some simple algorithms to determine animal movements with uncertainty based on light-level geolocators.


## Installing

`probGLS` now relies on `SGAT`, which can be installed from GitHub, using the devtools package. 

```R
devtools::install_github("SWotherspoon/SGAT")
```

`probGLS` itself can also be installed from GitHub, using the devtools package. 

```R
devtools::install_github("benjamin-merkel/probGLS")
```

If you don't have `devtools` installed already, install it first. 

```R
install.packages("devtools")
```

(probGLS or SGAT otherwise do not need devtools for normal use.)



## UPDATE Dec 2022

This package received a major update Dec 2022

- remove dependency on `raster` and `sp` and move to `sf`
- fix of a bug in `load_NOAA_OISST_V2`
- updated plotting functions to adhere with the updated `prob_algorithm` output
- updated manuals for each function