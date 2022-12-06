
#' Internal function to collapse lists into data.frame
#' 
#' Internal function to collapse lists into data.frame
#' @param x List that should be collapsed into a data.frame
#' @keywords internal


fun_list_to_dataframe  <- function(x) function(i) sapply(x, `[[`, i)

#' Internal function to calculate the proportion of time dry
#' 
#' Internal function to calculate the proportion of time between to twilight events a logger was registered as being dry
#' @param x point cloud data.frame
#' @param y act data.frame given to prob_algorithm
#' @param time_max End time point
#' @param wr sampling rate of conductivity switch in sec (e.g. MK15 & MK3006 sample every 3 sec) as given to prob_algorithm
#' @keywords internal


fun_time_dry <- function (x, y, time_max, wr) {
  if(!is.null(act)){
    wetdry            <- y$wetdry[y$dtime >= min(x$tFirst) & y$dtime <= time_max]
    wetdry_time_diff  <- abs(as.numeric(difftime(min(x$tFirst), time_max, units = 'secs')))
    sumact            <- (1 - sum(wetdry) * wr / wetdry_time_diff)
    if(sumact > 1) sumact <- 1
    if(sumact < 0) sumact <- 0
  } else {
    sumact <- 1
  }
  return(sumact)
}

#' Internal function to calculate the weight of each particle given the animals speed
#' 
#' Internal function to calculate the weight of each particle given the animals speed
#' @param x point cloud data.frame
#' @param spd speed.dry given to prob_algorithm
#' @param spw speed.wet given to prob_algorithm
#' @keywords internal


fun_weight_speed <- function(x, spd, spw) {
  
  m <- (spd[1] * x$prob.dry) + (spw[1] * (1 - x$prob.dry))
  s <- (spd[2] * x$prob.dry) + (spw[2] * (1 - x$prob.dry))
  v <- (spd[3] * x$prob.dry) + (spw[3] * (1 - x$prob.dry))
  
  w = dnorm(x$speed_ms, mean = m, sd = s)/max(dnorm(m, mean = m, sd = s), na.rm=T)
  
  w[x$speed_ms <= m] <- 1
  w[x$speed_ms >  v] <- 0
  w[x$speed_ms <  0] <- 0
  
  return(w)
}

#' Internal function to calculate the weight of each particle given SST and tag temperature
#' 
#' Internal function to calculate the weight of each particle given SST and tag temperature
#' @param x point cloud data.frame
#' @param y sst.sd given to prob_algorithm
#' @param z max.sst.diff given to prob_algorithm
#' @keywords internal


fun_weight_sst   <- function(x, y, z) {   
  w = dnorm(x$sst.diff, mean = 0, sd = y + x$sat.sst.err)/
      max(dnorm(0     , mean = 0, sd = y + x$sat.sst.err), na.rm=T)
  
  w[x$sst.diff >  z] <- 0 
  w[x$sst.diff < -z] <- 0 
  return(w)
}

#' Internal function to calculate the weight of each particle given sea ice concentration
#' 
#' Internal function to calculate the weight of each particle given sea ice concentration
#' @param x point cloud data.frame
#' @param y ice.conc.cutoff given to prob_algorithm
#' @keywords internal


fun_weight_ice   <- function(x, y) {   
  w <- rep(1, nrow(x))
  w[x$sat.ice > y] <- 0 
  return(w)
}


#' Internal function to extract on/off land mask from NOAA OISST V2
#' 
#' This internal function determines if points are on land using the NOAA OISST V2 land maks file.
#' @param FILE_NAME landmask file name with path
#' @param LONS Point longitudes
#' @param LATS Point latitudes
#' @keywords internal


load_landmask <- function(FILE_NAME, LONS, LATS){
  
  # open ncdf connection
  nc      <- nc_open(FILE_NAME)
  nclon   <- ncvar_get(nc, "lon")
  nclat   <- ncvar_get(nc, "lat")
  
  LONS[LONS < 0] <- LONS[LONS < 0] + 360
  
  lonindx <- apply(matrix(LONS), 1, function(x) which.min(abs(x - nclon)))
  latindx <- apply(matrix(LATS), 1, function(x) which.min(abs(x - nclat)))
  
  ncdat   <- mapply(function(x, y) {
    ncvar_get(nc, 
              varid = 'lsmask',
              start = c(x, y, 1),
              count = c(1, 1, 1))
  }, x = lonindx, y = latindx, SIMPLIFY = F)
  
  ncdat   <- unlist(ncdat)
  
  # close ncdf connection
  nc_close(nc)
  
  return(ncdat)
}
