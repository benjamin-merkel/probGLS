
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


#' Calculate MCPs
#' 
#' This function calculates minimum convex polygons (MCP) without calculating area and without minimum required points. IT is independent of adehabitatHR.
#' @param xy data in x dimension to be plotted
#' @param percent Percentage considered to determine MCP
#' @param prob1 lower quantile considered
#' @param prob2 upper quantile considered

calc_mcp <- function (xy, percent = 100) {
  
  id      <- xy[[1]]
  id      <- factor(id)
  xy      <- as.data.frame(st_coordinates(xy))
  r       <- split(xy, id)
  est.cdg <- function(xy) apply(xy, 2, mean)
  cdg     <- lapply(r, est.cdg)
  levid   <- levels(id)
  res     <- lapply(1:length(r), function(i) {
    k        <- levid[i]
    df.t     <- r[[levid[i]]]
    cdg.t    <- cdg[[levid[i]]]
    dist.cdg <- function(xyt) {
      d <- sqrt(((xyt[1] - cdg.t[1])^2) + ((xyt[2] - cdg.t[2])^2))
      return(d)
    }
    di       <- apply(df.t, 1, dist.cdg)
    key      <- c(1:length(di))
    acons    <- key[di <= quantile(di, percent/100)]
    xy.t     <- df.t[acons, ]
    coords.t <- chull(xy.t[, 1], xy.t[, 2])
    xy.bord  <- xy.t[coords.t, ]
    xy.bord  <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    so       <- st_as_sf(xy.bord, coords = c("X", "Y"), crs = 4326)
    so       <- st_combine(so)
    so       <- st_cast(so, "POLYGON")
    so       <- st_sf(so)
    return(so)
  })
  df <- do.call(rbind, res)
  df$id <- levels(id)
  
  return(df)
}

