#' load NOAA OISST V2
#' 
#' This function extracts information from NOAA OISST V2 NetCDF files available at
#' ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/ directory
#' @param FILE_NAME full path to NetCDF data file
#' @param LONS Point longitudes
#' @param LATS Point latitudes
#' @param DATE date considered
#' @param extract.value which data to extract: "sst" - SST, "err" - SST error, "icec" - sea ice concentration
#' @return A vector of values for each lon lat pair given the extract.value
#' @return  
#' @return NetCDF files should be downloaded from the links on:
#' @return http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
#' @return In addition to the temperature data files, also download a copy of the landmask file lsmask.oisst.v2.nc from the same page. Inside the NetCDF files, data are available on a 0.25 degree latitude x 0.25 degree longitude global grid (720x1440 cells) From -89.875N to 89.875N, 0.125E to 359.875E. Locations are at the CENTER of a grid cell.
#' @export


load_NOAA_OISST_V2 <- function(FILE_NAME, LONS, LATS, DATE, extract.value = 'sst'){
  
  # open ncdf connection
  nc      <- nc_open(FILE_NAME)
  
  nclon   <- ncvar_get(nc, "lon")
  nclat   <- ncvar_get(nc, "lat")
  ncdates <- as.Date(ncvar_get(nc, "time"), origin = '1800-1-1') #available time points in nc
  
  LONS[LONS < 0] <- LONS[LONS < 0] + 360
  
  lonindx <- apply(matrix(LONS), 1, function(x) which.min(abs(x - nclon)))
  latindx <- apply(matrix(LATS), 1, function(x) which.min(abs(x - nclat)))
  dateindx<- which.min(abs(DATE - ncdates))
  
  ncdat   <- mapply(function(x, y, z) {
              ncvar_get(nc, 
                        varid = extract.value,
                        start = c(x, y, z),
                        count = c(1, 1, 1))
  }, x = lonindx, y = latindx, z = dateindx, SIMPLIFY = F)
  
  ncdat   <- unlist(ncdat)
  
  # close ncdf connection
  nc_close(nc)
  
  # If there are missing data in the NetCDF, they should appear as 32767.
  # Replace that value with NA if it occurs anywhere.
  ncdat = ifelse(ncdat == 32767, NA, ncdat)
  
  return(ncdat)
  
}
