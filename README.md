# probGLS

The package probabilistic algorithm for geolocation data (probGLS) provides some simple algorithms to determine animal movements with uncertainty based on light-level geolocators.


## Installing

`probGLS` now relies on `SGAT` (https://github.com/SWotherspoon/SGAT) and `GeoLight` (https://github.com/slisovski/GeoLight). These can be installed from GitHub, using the devtools package. 

```R
devtools::install_github("SWotherspoon/SGAT")
devtools::install_github("SLisovski/GeoLight")
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

This package received a major update Dec 2022.

- removed dependency on `raster` and `sp` (excluding `spDists` used when distance method chosen is "ellipsoid" bc it is faster than `st_distance`) and move to `sf`
- fixed a bug in `load_NOAA_OISST_V2` and streamlined the function
- updated plotting functions to adhere with the updated `prob_algorithm` output
- updated manuals for each function