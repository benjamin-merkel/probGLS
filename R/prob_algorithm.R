#' probabilistic algorithm for geolocation data
#' 
#' prob_algorithm is a simple probabilistic algorithm to be used with geolocation data.
#' @param particle.number number of particles for each location cloud used in the model
#' @param iteration.number number of iterations
#' @param loess.quartile quartiles for loessFilter (GeoLight), if NULL loess filter is not used
#' @param tagging.location tagging location longitude and latitude
#' @param tagging.date deployment data as POSIXct or Date object
#' @param retrieval.date  retrieval date as POSIXct or Date object
#' @param sunrise.sd output vector from twilight_error_estimation
#' @param sunset.sd output vector from twilight_error_estimation 
#' @param tol 
#' @param range.solar min and max of solar angle range in degree
#' @param speed.dry optimal speed, speed standard deviation and max speed allowed if logger is dry in m/s
#' @param speed.wet optimal speed, speed standard deviation and max speed allowed if logger is wet in m/s
#' @param distance.method  
#' @param sst.sd SST standard deviation in degree C 
#' @param max.sst.diff max difference in SST allowed in degree C
#' @param ice.conc.cutoff max percentage of sea ice in which the animal is believed to be
#' @param boundary.box min lon, max lon, min lat and max lat of extrem boundary where you expect an animal to be
#' @param east.west.comp if T apply biotrack east west movement compensation (Biotrack manual v11 page 31pp.)
#' @param land.mask if T animal is only using ocean areas, if F animal is only using land areas, if NULL no land mask used
#' @param med.sea if T classifiy mediterranean sea as land   
#' @param black.sea if T classifiy black sea as land      
#' @param baltic.sea if T classifiy baltic sea as land   
#' @param caspian.sea if T classifiy caspian sea as land   
#' @param sensor data.frame with daily SST data deduced from tag temperature readings (sst_deduction ouput), NULLif no SST data is available (SST will not be used)
#' @param trn data.frame containing twilights and at least tFirst, tSecond and type (same as computed by trn_to_dataframe, ipe_to_dataframe or lotek_to_dataframe)
#' @param act data.frame containing wet dry data (e.g. .act file from Biotrack loggers or .deg file from migrate tech loggers), NULL if no wetdry data is available (algorithm will assume that the logger was always dry)
#' @param wetdry.resolution sampling rate of conductivity switch in sec (e.g. MK15 & MK3006 sample every 3 sec)
#' @param backward run algorithm from end to start
#' @param NOAA.OI.location directory location of NOAA OI V2 NCDF files as well as land mask file 'lsmask.oisst.v2.nc' (downloadable from http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html)
#' @return A list with: [1] all positions, [2] geographic median positions, [3] all possible particles, [4] input parameters, [5] model run time. List items 1 to 3 are returned as SpatialPointsDataframe.
#' @details Many weighting parameters can be used. Some others (which are not yet implemented) are: surface air temperature, air pressure, water salinity, and topography/ bathymetry.
#' @import ncdf4
#' @ImportFrom sp spDists
#' @import GeoLight
#' @import SGAT
#' @import sf
#' @examples
#'######################################
#'# example black browed albatross track 
#'######################################
#' 
#'# define start and end datetimes ----
#'start <- as.POSIXct("2014-12-13 17:55", tz="UTC")
#'end   <- as.POSIXct("2014-12-22 08:55", tz="UTC")
#'
#'# light data ----
#'trn           <- twilightCalc(BBA_lux$dtime, BBA_lux$lig, ask = FALSE, LightThreshold = 2, maxLight = 5)
#'
#'# sst data ----
#'sen           <- sst_deduction(datetime = BBA_sst$dtime, temp = BBA_sst$temp, temp.range = c(-2,30))
#'
#'# wet dry data ----
#'act           <- BBA_deg[BBA_deg$wet.dry=="wet",]
#'act$wetdry    <- act$duration
#'
#'# twilight error distribution estimation ----
#'tw            <- twilight_error_estimation()
#'
#'# download environmental data ----
#'
#'# download yearly NetCDF files for (replace YEAR with appropriate number): 
#'# daily mean SST                   -> 'sst.day.mean.YEAR.v2.nc'
#'# daily SST error                  -> 'sst.day.err.YEAR.v2.nc'
#'# daily mean sea ice concentration -> 'icec.day.mean.YEAR.v2.nc'
#'
#'# from:
#'# https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
#'# and place all into the same folder
#'
#'# Also, download the land mask file: 'lsmask.oisst.v2.nc' from the same directory 
#'# and place it in the same folder as all the other NetCDF files
#'
#'# run algorithm ----
#'pr   <- prob_algorithm(trn                         = trn, 
#'                       sensor                      = sen[sen$SST.remove==F,], 
#'                       act                         = act, 
#'                       tagging.date                = min(trn$tFirst), 
#'                       retrieval.date              = max(trn$tSecond), 
#'                       loess.quartile              = NULL, 
#'                       tagging.location            = c(-36.816,-54.316), 
#'                       particle.number             = 2000, 
#'                       iteration.number            = 100,
#'                       sunrise.sd                  = tw,
#'                       sunset.sd                   = tw,
#'                       range.solar                 = c(-7,-1),
#'                       boundary.box                = c(-120,40,-90,0),
#'                       speed.dry                   = c(12,6,45),
#'                       speed.wet                   = c(1,1.3,5), 
#'                       sst.sd                      = 0.5,       
#'                       max.sst.diff                = 3,          
#'                       east.west.comp              = T,
#'                       land.mask                   = T, 
#'                       ice.conc.cutoff             = 1, 
#'                       tol                         = 0.08,
#'                       wetdry.resolution           = 1,
#'                       distance.method             = "ellipsoid"
#'                       NOAA.OI.location            = "folder with environmental data and land mask")
#'
#'# plot lat, lon, SST vs time ----
#'plot_timeline(pr, solar.angle = mean(pr[[2]]$median.solar.angle))
#'
#'# plot lon vs lat map ----
#'plot_map(pr, legend.position = "topright")
#' @export


prob_algorithm <- function(
    particle.number             = 100,
    iteration.number            = 10,
    loess.quartile              = NULL,
    trn                         = trn,
    sensor                      = sen,
    act                         = act,
    tagging.date                = start, 
    retrieval.date              = end,
    tol                         = 0.08,
    tagging.location            = c(-36.816,-54.316),
    boundary.box                = c(-180,180,-90,90),
    sunrise.sd                  = tw,
    sunset.sd                   = tw,
    range.solar                 = c(-7,-1),
    speed.wet                   = c(20,0.2,25),
    speed.dry                   = c(20,0.2,25),
    sst.sd                      = 0.5,
    max.sst.diff                = 3,
    # days.around.spring.equinox  = c(10,10),
    # days.around.fall.equinox    = c(10,10),
    ice.conc.cutoff             = 1,
    wetdry.resolution           = 1,
    east.west.comp              = T,  
    land.mask                   = T,
    med.sea                     = T,
    black.sea                   = T,
    baltic.sea                  = T,
    caspian.sea                 = T,
    backward                    = F,
    distance.method             = "ellipsoid", # c("spherical", "ellipsoid") spherical is slower but uses s2, ellipsoid is much faster and uses great circle on ellipsoid
    NOAA.OI.location            = 'E:/environmental data/SST/NOAA OI SST V2'){
  
  start.time <- Sys.time()
  
  # appease R CMD check
  tFirst <- tSecond <- type <- dtime <- doy <- jday <- year <- month <- NULL
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  
  if(is.null(sensor)) sst.used=F else sst.used=T
  model.input <- data.frame(parameter=c('particle.number','iteration.number','loess.quartile','tagging.location',
                                        'tagging.date','retrieval.date','sunrise.sd','sunset.sd','range.solar','speed.wet',
                                        'speed.dry','sst.sd','max.sst.diff','days.around.spring.equinox',
                                        'days.around.fall.equinox','ice.conc.cutoff','boundary.box','med.sea','black.sea',
                                        'baltic.sea','caspian.sea','east.west.comp','wetdry.resolution','NOAA.OI.location','backward','sensor.data'),
                            chosen=c(paste(particle.number,collapse=" "),paste(iteration.number,collapse=" "),paste(loess.quartile,collapse=" "),paste(tagging.location,collapse=" "),
                                     paste(tagging.date,collapse=" "),paste(retrieval.date,collapse=" "),paste(sunrise.sd,collapse=" "),paste(sunset.sd,collapse=" "),paste(range.solar,collapse=" "),paste(speed.wet,collapse=" "),
                                     paste(speed.dry,collapse=" "),paste(sst.sd,collapse=" "),paste(max.sst.diff,collapse=" "),paste(days.around.spring.equinox,collapse=" "),
                                     paste(days.around.fall.equinox,collapse=" "),paste(ice.conc.cutoff,collapse=" "),paste(boundary.box,collapse=" "),paste(med.sea,collapse=" "),paste(black.sea,collapse=" "),
                                     paste(baltic.sea,collapse=" "),paste(caspian.sea,collapse=" "),paste(east.west.comp,collapse=" "),paste(wetdry.resolution,collapse=" "),paste(NOAA.OI.location,collapse=" "),
                                     paste(backward,collapse=" "),sst.used))
  
  # test if land mask file is available if chosen ----
  if(!is.null(land.mask)){
    landmask.location <- list.files(path=NOAA.OI.location,pattern="lsmask.oisst.v2.nc",recursive=T)
    if(length(landmask.location)==0){
      # stop(paste('no land mask file found in folder',NOAA.OI.location,sep=' '),call.=F)
      cat('\r','no land mask file found - file will be downloaded from noaa.gov')
      download.file('https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/lsmask.oisst.v2.nc', 'data/lsmask.oisst.v2.nc', mode="wb")
      landmask.location <- paste0(getwd(),"/data")
    } else {
      landmask.location <- paste(NOAA.OI.location,landmask.location,sep='/')[1]
    }
  }
  
  # add date time julian doy etc -----
  trn$dtime     <- trn$tFirst + as.numeric(difftime(trn$tSecond,trn$tFirst,units='sec'))/2 # mean between tFirst and tSecond
  trn$doy       <- as.numeric(strftime(trn$dtime, format = "%j"))
  trn$month     <- as.numeric(strftime(trn$dtime, format = "%m"))
  trn$year      <- as.numeric(strftime(trn$dtime, format = "%Y"))
  trn$jday      <- as.numeric(julian(trn$dtime))
  
  # remove outside known data -----
  trn    <- trn[trn$tFirst >= as.POSIXct(tagging.date) & trn$tSecond <= as.POSIXct(retrieval.date),]
  trn    <- trn[!is.na(trn$tFirst) & !is.na(trn$tSecond),]
  if(nrow(trn)==0)  stop('no data points in trn file between selected tagging and retrieval date',call.=F)
  
  # remove all twilight outliers ----
  if(!is.null(loess.quartile)){
    trn$loess  <- loessFilter(trn, plot = T, k = loess.quartile)
    trn        <- trn[trn$loess == T,]
  }

  # west east movement compensation----
  if(east.west.comp==T){
    timedate   <- trn$dtime
    temp.lon   <- as.data.frame(coord(trn, degElevation = -6, note = F, method = 'NOAA'))[,1]
    
    lon2       <- -temp.lon
    lon1       <- -c(NA, temp.lon[1:(length(temp.lon) - 1)])
    timedate2  <- as.numeric(timedate)
    timedate1  <- as.numeric(c(NA, timedate[1:(length(timedate) - 1)]))
    length     <- abs(as.numeric(difftime(trn$tSecond,trn$tFirst,units = "sec")))
    m          <- ((lon2 - lon1) / (timedate2 - timedate1)) * (length / 15) /2
    ewc        <- data.frame(timedate, lon2, lon1, timedate2, timedate1, length, m)
    
    # add/substract m to sunsets/sunrises in hours
    ewc$tFirst.corrected  <- trn$tFirst  + ewc$m * 3600
    ewc$tSecond.corrected <- trn$tSecond - ewc$m * 3600
    ewc$tFirst.corrected [is.na(ewc$tFirst.corrected)]  <- trn$tFirst [is.na(ewc$tFirst.corrected)]
    ewc$tSecond.corrected[is.na(ewc$tSecond.corrected)] <- trn$tSecond[is.na(ewc$tSecond.corrected)]
    
    trn$tFirst  <- ewc$tFirst.corrected
    trn$tSecond <- ewc$tSecond.corrected
  }
  
  # subset data to columns of interest-----
  trn     <- data.frame(subset(trn, select=c(tFirst, tSecond, type, dtime, doy, jday, year, month)))
  
  # duplicate data frame number of particle times----
  trn     <- data.frame(mapply(rep, trn, particle.number))
  
  # make sure date time columns are in the right format-----
  trn$'tFirst'  <- as.POSIXct(as.numeric(as.character(trn$'tFirst')),  origin="1970-01-01", tz="UTC")
  trn$'tSecond' <- as.POSIXct(as.numeric(as.character(trn$'tSecond')), origin="1970-01-01", tz="UTC")
  trn$'dtime'   <- as.POSIXct(as.numeric(as.character(trn$'dtime')),   origin="1970-01-01", tz="UTC")
  
  # assign unique step to each particle cloud -----
  trn$step      <- as.numeric(as.factor(paste(as.numeric(julian(trn$tFirst)),as.numeric(julian(trn$tSecond)),trn$type,sep="-") )) 
  
  # add sun elevation angle -----
  solar.angle.steps <- seq(range.solar[1], range.solar[2], 0.01)
  trn$solar.angle   <- sample(solar.angle.steps, size = nrow(trn), replace = T)
  
  # vary tFirst and tSecond ----
  trn$tFirst.err [trn$type==1] <-  60 * (rlnorm(length(trn[trn$type==1,1]), meanlog = sunrise.sd[1], sdlog = sunrise.sd[2]) + sunrise.sd[3])
  trn$tSecond.err[trn$type==1] <- -60 * (rlnorm(length(trn[trn$type==1,1]), meanlog = sunset.sd[1],  sdlog = sunset.sd[2])  + sunset.sd[3]) 
  trn$tFirst.err [trn$type==2] <- -60 * (rlnorm(length(trn[trn$type==2,1]), meanlog = sunset.sd[1],  sdlog = sunset.sd[2])  + sunset.sd[3]) 
  trn$tSecond.err[trn$type==2] <-  60 * (rlnorm(length(trn[trn$type==2,1]), meanlog = sunrise.sd[1], sdlog = sunrise.sd[2]) + sunrise.sd[3]) 
  trn$tFirst     <- trn$tFirst  + trn$tFirst.err
  trn$tSecond    <- trn$tSecond + trn$tSecond.err
  
  # calculate coordinates out of twilight times and sun elevation angle----
  pos <- thresholdEstimate(trise  = ifelse(trn$type == 1, trn$tFirst, trn$tSecond),
                           tset   = ifelse(trn$type == 1, trn$tSecond, trn$tFirst),
                           zenith = 90 - trn$solar.angle,
                           tol    = tol)
  trn$lon <- pos[,1]
  trn$lat <- pos[,2]
  
  # remove sun elevation angle during equinox----
  # trn$solar.angle[trn$doy %in% c(spring.equinox,fall.equinox)] <- NA
  trn$solar.angle[is.na(trn$lat)] <- NA
  
  # define equinox periods----
  # doy 79  = 20 March
  # doy 265 = 22 September
  # spring.equinox <- c((79  - days.around.spring.equinox[1]):(79  + days.around.spring.equinox[2]))
  # fall.equinox   <- c((265 - days.around.fall.equinox[1])  :(265 + days.around.fall.equinox[2]))
  
  # assign random latitudes to equinox periods----
  # trn$lat[is.na(trn$lat) & trn$doy %in% c(spring.equinox,fall.equinox)] <- sample(seq(boundary.box[3],boundary.box[4],by=0.0001),size=length(trn$lat[is.na(trn$lat) & trn$doy %in% c(spring.equinox,fall.equinox)]),replace=T)
  trn$lat[is.na(trn$lat)] <- sample(seq(boundary.box[3],boundary.box[4],by=0.0001),size=length(trn$lat[is.na(trn$lat)]), replace=T)
  
  
  # remove all positions outside boundaries-----
  if(boundary.box[1] > boundary.box[2]){
    trn    <- trn[(trn$lon>boundary.box[1] | trn$lon<boundary.box[2]) & trn$lat>boundary.box[3] & trn$lat<boundary.box[4],]
  } else {
    trn    <- trn[ trn$lon>boundary.box[1] & trn$lon<boundary.box[2]  & trn$lat>boundary.box[3] & trn$lat<boundary.box[4],]
  }
  
  # remove all point clouds smaller then 1/5 of particles ----
  trn              <- trn[!is.na(trn$lat),]
  jt               <- data.frame(table(trn$step))
  trn              <- trn[trn$step %in% as.numeric(as.character(jt$Var1[jt$Freq >= c(particle.number/5)])),]
  
  # remove point clouds with maximum lower then min.lat ----
  rm.lat           <- data.frame(max.lat = tapply(trn$lat, trn$step, max))
  rm.lat$step      <- rownames(rm.lat)
  
  if(nrow(rm.lat)==0){
    stop(paste("No data points inside boundary box. increase boundary box"),call.=F)  
  }
  
  rm.lat$rm        <- 1
  rm.lat$rm[rm.lat$max.lat < boundary.box[3]] <- 0
  trn              <- trn[trn$step %in% rm.lat$step[rm.lat$rm==1],]
  
  # remove point clouds with minimum higher then max.lat----
  rm.lat           <- data.frame(min.lat=tapply(trn$lat, trn$step, min))
  rm.lat$step      <- rownames(rm.lat)
  
  if(nrow(rm.lat)==0){
    stop(paste("No data points inside boundary box. increase boundary box"),call.=F)  
  }
  
  rm.lat$rm        <- 1
  rm.lat$rm[rm.lat$min.lat > boundary.box[4]] <- 0
  trn              <- trn[trn$step %in% rm.lat$step[rm.lat$rm==1],]
  
  
  # remove everything on land-----
  if(!is.null(land.mask)) trn$landmask <- load_landmask(FILE_NAME = landmask.location, LONS = trn$lon, LATS = trn$lat)
  # if(!is.null(land.mask))  landms       <- rotate(raster(landmask.location))
  # if(!is.null(land.mask))  trn$landmask <- extract(landms, trn)
  
  
  # remove baltic sea ---- 
  if(baltic.sea==T)  trn$landmask[trn$lon>14 & trn$lon<33.5 & trn$lat>51.4 & trn$lat<66.2] <- 0
  
  # remove meditereanian sea -----
  if(med.sea==T) {
    trn$landmask[trn$lon>=0  & trn$lon<=27  & trn$lat>30 & trn$lat<48] <- 0
    trn$landmask[trn$lon>=27 & trn$lon< 40  & trn$lat>30 & trn$lat<40] <- 0
    trn$landmask[trn$lon>355 & trn$lon<=360 & trn$lat>30 & trn$lat<42] <- 0
  }
  
  # remove black sea ---- 
  if(black.sea==T)   trn$landmask[trn$lon>27 & trn$lon<45 & trn$lat>40 & trn$lat<48] <- 0
  
  # remove caspian sea ---- 
  if(caspian.sea==T) trn$landmask[trn$lon>45 & trn$lon<62 & trn$lat>35 & trn$lat<48] <- 0
  
  
  # remove all point clouds with less than 10 % points outside land ----
  if(!is.null(land.mask)){
    if(land.mask==T) jt <- trn[trn$landmask==1,]
    if(land.mask==F) jt <- trn[trn$landmask==0,]
    
    jt               <- data.frame(table(jt$step))
    
    trn              <- trn[trn$step %in% jt$Var1[jt$Freq >= particle.number * 0.1],]
    
    trn$weight_land <- 0
    if(land.mask==T) trn$weight_land[trn$landmask==1] <- 1
    if(land.mask==F) trn$weight_land[trn$landmask==0] <- 1
  }
  
  
  # define date and time as average between tFirst and tSecond-----
  trn$dtime        <- trn$tFirst + as.numeric(difftime(trn$tSecond,trn$tFirst,units='sec'))/2
  trn$jday2        <- floor(trn$jday)
  trn$date         <- as.Date(trn$dtime)
  trn              <- trn[order(trn$dtime),]
  
  # transform to SF object -----
  trn              <- st_as_sf(trn, coords = c("lon", "lat"), crs = 4326)
  trn              <- cbind(trn,lon = st_coordinates(trn)[,1], lat = st_coordinates(trn)[,2])
  
  # create data frame of all particles computed ----
  all.particles    <- trn
  
  
  # create spatial tagging.location object-----
  tag.loc       <- data.frame(lon = tagging.location[1], lat = tagging.location[2])
  tag.loc       <- st_as_sf(tag.loc, coords = c("lon", "lat"), crs = 4326)
  tag.loc       <- cbind(tag.loc, lon = st_coordinates(tag.loc)[,1], lat = st_coordinates(tag.loc)[,2], 
                         dtime = as.POSIXct(tagging.date, tz="UTC"), date = as.Date(tagging.date), 
                         doy = NA, jday = NA, jday2 = NA, year = NA, month = NA, type = NA, prev.dtime = NA,
                         tFirst = as.POSIXct(tagging.date,tz="UTC"), tSecond = as.POSIXct(tagging.date,tz="UTC"),
                         tFirst.err = NA, tSecond.err = NA, solar.angle = NA, step = NA, iteration = NA, landmask = NA,
                         dist_m = NA, prob.dry = NA, speed_ms = NA, time.diff = NA, sat.ice = NA, sat.sst = NA,
                         sat.sst.err = NA, tag.sst = NA, sst.diff = NA, weight_land = 1, weight_ice = NA, weight_sst = NA, 
                         weight_speed = NA, rel_weight = 1)
  
  
  # create current location list and start with tagging.location and time of tagging-----
  loc             <- tag.loc
  if(backward==F) loc$dtime       <- as.POSIXct(tagging.date)
  if(backward==T) loc$dtime       <- as.POSIXct(retrieval.date)
  loc$jday        <- as.numeric(julian(loc$dtime))
  
  loc.list        <- vector("list", length = iteration.number)
  loc.list        <- lapply(loc.list, function (x) st_as_sf(loc))
  
  # loop through each step ----
  iter = 0
  
  if(backward==F) steps <- sort(unique(trn$step))
  if(backward==T) steps <- sort(unique(trn$step),decreasing=T)
  
  progress_bar = txtProgressBar(min=0, max=length(steps), style = 3, char='=')
  
  
  for(ts in steps){
    step.start <- Sys.time()
    
    if(length(trn$dtime[trn$step == ts]) > 0){
      
      # select the right time step----
      trn.step      <- trn[trn$step == ts,]
      
      # calculate bearing and distance to previous locations----
      if(distance.method == "spherical"){
        suppressMessages(sf_use_s2(T))
        gdist <- lapply(loc.list, FUN = function(x) as.numeric(st_distance(trn.step, x)))  # in meters
      }
      
      if(distance.method == "ellipsoid"){
        gdist <- lapply(loc.list, FUN=function(x) spDists(matrix(c(trn.step$lon, trn.step$lat), ncol = 2), 
                                                          matrix(c(x$lon, x$lat), ncol = 2), longlat = T) * 1000)
      }
      
      gprob.dry    <- lapply(loc.list, fun_time_dry, y = act, time_max = max(trn.step$tSecond), wr = wetdry.resolution)  
      
      # create list with one data frame for each iteration-----
      gr2 <- vector("list", length = iteration.number)
      gr2 <- lapply(gr2, function (x) trn.step)
      it2 <- as.list(1:iteration.number)
      
      # calculate speed and time difference-----
      prev.dtime <- lapply(loc.list, function(x) x$dtime)
      gr2        <- mapply(cbind, gr2, prev.dtime = prev.dtime, iteration = it2, SIMPLIFY = FALSE)
      time.diff  <- lapply(gr2, function(x) as.numeric(difftime(x$dtime, x$prev.dtime, units="secs")))
      gr2        <- mapply(cbind, gr2, dist_m = gdist, prob.dry = gprob.dry, time.diff = time.diff, SIMPLIFY = FALSE)
      speed_ms   <- lapply(gr2,function(x) x$dist_m / x$time.diff)
      gr2        <- mapply(cbind, gr2, speed_ms = speed_ms, SIMPLIFY = FALSE)
      
      # load needed SST, SST error and ice data-----
      ls <- list.files(NOAA.OI.location)
      fname.sst <- paste0(NOAA.OI.location,'/',ls[grep(paste0('sst.day.mean.'  ,trn.step$year[1]),ls)])[1]
      fname.err <- paste0(NOAA.OI.location,'/',ls[grep(paste0('sst.day.err.'   ,trn.step$year[1]),ls)])[1]
      fname.ice <- paste0(NOAA.OI.location,'/',ls[grep(paste0('icec.day.mean.' ,trn.step$year[1]),ls)])[1]
      
      trn.step$sat.sst     <- load_NOAA_OISST_V2(FILE_NAME = fname.sst, LONS = trn.step$lon, LATS = trn.step$lat, DATE = trn.step$date[1], extract.value = 'sst')
      trn.step$sat.ice     <- load_NOAA_OISST_V2(FILE_NAME = fname.ice, LONS = trn.step$lon, LATS = trn.step$lat, DATE = trn.step$date[1], extract.value = 'icec')
      trn.step$sat.sst.err <- load_NOAA_OISST_V2(FILE_NAME = fname.err, LONS = trn.step$lon, LATS = trn.step$lat, DATE = trn.step$date[1], extract.value = 'err')
      trn.step$sat.ice[is.na(trn.step$sat.ice)] <- 0
      
      # remove sst values in pixels with more than ice.conc.cutoff----
      # except for 2012-8-11 to 2012-8-16 as ice data is screwed up in this period
      if(as.Date(trn.step$date[1])<as.Date("2012-08-11") | 
         as.Date(trn.step$date[1])>as.Date("2012-08-16"))  trn.step$sat.sst[trn.step$sat.ice > ice.conc.cutoff] <-NA
      
      
      # set tag temp reading to NA if no sensor data was available 
      if(!is.null(sensor)) sensor$jday  <- as.numeric(julian(sensor$date))
      if(length(sensor$SST[sensor$jday == trn.step$jday2[1]]) > 0)  trn.step$tag.sst <- sensor$SST[sensor$jday == trn.step$jday2[1]]
      if(length(sensor$SST[sensor$jday == trn.step$jday2[1]]) ==0)  trn.step$tag.sst <- NA
      
      
      trn.step$sst.diff <- trn.step$sat.sst - trn.step$tag.sst
      
      gr2   <- lapply(gr2,function(x) cbind(x,
                                            sat.sst     = trn.step$sat.sst,
                                            tag.sst     = trn.step$tag.sst,
                                            sst.diff    = trn.step$sst.diff,
                                            sat.sst.err = trn.step$sat.sst.err,
                                            sat.ice     = trn.step$sat.ice))
      if(!is.null(land.mask)) if(land.mask==T) gr2   <- lapply(gr2,function(x) x[!is.na(x$sat.sst),])
      
      
      # weigh each particle according to speed, ice and SST
      wspeed <- lapply(gr2, fun_weight_speed, spd = speed.dry, spw = speed.wet)  
      wsst   <- lapply(gr2, fun_weight_sst, y = sst.sd, z = max.sst.diff)  
      wice   <- lapply(gr2, fun_weight_ice, y = ice.conc.cutoff)  
      
      
      # if no tag.SST available only use speed weighing----
      if(!is.na(trn.step$tag.sst[1])) wrel  <- mapply(function(x, y, z) {return(x * y * z)}, x = wspeed, y = wsst, z = wice, SIMPLIFY = F)
      if( is.na(trn.step$tag.sst[1])) wrel  <- mapply(function(x, z) {return(x * z)}, x = wspeed, z = wice, SIMPLIFY = F)
      
      gr2    <- mapply(cbind, gr2, weight_speed = wspeed, weight_sst = wsst, weight_ice = wice, rel_weight = wrel, SIMPLIFY = FALSE)
      if(!is.null(land.mask)) lapply(gr2,function(x) {x$rel_weight <- x$rel_weight * x$weight_land})
      
      
      # choose a specific random point based on weighing of particles ----
      # if all weights are 0 and/or NA jump this position
      chosen.point   <- vector("list", length = iteration.number)
      new.loc        <- vector("list", length = iteration.number)
      
      chosen.point   <- lapply(gr2, function(x) {
        if(max(x$rel_weight, na.rm=T) > 0) {
          return(sample(1:nrow(x), size = 1, prob = x$rel_weight))
        } else {
          return(0)
        }
      })
      
      new.loc        <- mapply(function(x, y, z) {
        if(y %in% 1:particle.number) {
          return(x[y,])
        } else {
          x <- z
          x$rel_weight <- 2
          return(x)
        }
      }, x = gr2, y = chosen.point, z = loc.list, SIMPLIFY = F)
      
      
      # save loop step data in object ----
      loc.list <- new.loc
      
      iter = iter + 1
      new.loc.df <- lapply(new.loc, function(x) data.frame(st_drop_geometry(x)))
      new.loc.df <- data.frame(Map(fun_list_to_dataframe(new.loc.df), names(new.loc.df[[1]])))
      
      if(iter==1) new.loc.df2 <- new.loc.df else new.loc.df2 <- rbind(new.loc.df2, subset(new.loc.df, select = names(new.loc.df2)))
      
    }
    
    step.end  <- Sys.time()
    step.time <- step.end - step.start
    
    # cat('\r',paste(as.Date(trn.step$dtime)[1],'  -  ',iter," of ",length(unique(trn$step)),' steps      ',sep=""))
    setTxtProgressBar(progress_bar, value = ts)
  }
  
  
  # remove all empty steps
  new.loc.df2        <- new.loc.df2[new.loc.df2$rel_weight <= 1,]
  
  #remove NAs in lat and lon
  new.loc.df2$lon <- as.numeric(new.loc.df2$lon)
  new.loc.df2$lat <- as.numeric(new.loc.df2$lat)
  new.loc.df2     <- new.loc.df2[!is.na(new.loc.df2$lon),]
  new.loc.df2     <- new.loc.df2[!is.na(new.loc.df2$lat),]
  
  new.loc.sf              <- st_as_sf(new.loc.df2, coords = c("lon", "lat"), crs = 4326)
  new.loc.sf$dtime        <- as.POSIXct(new.loc.sf$dtime  ,origin="1970-01-01",tz="UTC")
  new.loc.sf$tFirst       <- as.POSIXct(new.loc.sf$tFirst ,origin="1970-01-01",tz="UTC")
  new.loc.sf$tSecond      <- as.POSIXct(new.loc.sf$tSecond,origin="1970-01-01",tz="UTC")
  new.loc.sf$jday2        <- as.numeric(julian(new.loc.sf$dtime))
  new.loc.sf              <- new.loc.sf[order(new.loc.sf$step),]
  
  # split points by "step" id
  sfsplit  <- split(new.loc.sf, new.loc.sf$step)
  
  # cast them to multipoints and combine
  sfcomb   <- do.call(c, Map(st_combine, sfsplit))
  
  # calculate geographic median for each particle cloud----
  median.loc.sf                    <- st_as_sf(st_point_on_surface(sfcomb))
  median.loc.sf$dtime              <- as.POSIXct(tapply(new.loc.df2$dtime, new.loc.df2$step, mean), origin="1970-01-01", tz="UTC")
  median.loc.sf$tFirst             <- as.POSIXct(tapply(new.loc.df2$tFirst, new.loc.df2$step, median), origin="1970-01-01", tz="UTC")
  median.loc.sf$tSecond            <- as.POSIXct(tapply(new.loc.df2$tSecond, new.loc.df2$step, median), origin="1970-01-01", tz="UTC")
  median.loc.sf$type               <- tapply(new.loc.df2$type, new.loc.df2$step, median)
  median.loc.sf$tag.sst            <- tapply(new.loc.df2$tag.sst, new.loc.df2$step, mean)
  median.loc.sf$mean.sat.sst       <- tapply(new.loc.df2$sat.sst, new.loc.df2$step, mean)
  median.loc.sf$median.solar.angle <- tapply(new.loc.df2$solar.angle, new.loc.df2$step, median)
  median.loc.sf$wmean.solar.angle  <- sapply(split(new.loc.df2, new.loc.df2$step), function(x) weighted.mean(x$solar.angle, x$rel_weight))
  median.loc.sf$mean.weight_sst    <- tapply(new.loc.df2$weight_sst, new.loc.df2$step, mean)
  median.loc.sf$mean.weight_speed  <- tapply(new.loc.df2$weight_speed, new.loc.df2$step, mean)
  median.loc.sf$mean.rel_weight    <- tapply(new.loc.df2$rel_weight, new.loc.df2$step, mean)
  
  sapply(split(new.loc.df2, new.loc.df2$step), function(x) weighted.mean(x$solar.angle, x$rel_weight))
    
  end.time   <- Sys.time()
  time.taken <- abs(difftime(end.time,start.time,units="mins"))
  
  cat('\r',paste('algorithm run time:',round(as.numeric(time.taken),1),'min      ',sep=" "))
  
  list.all            <- list(new.loc.sf, median.loc.sf, all.particles, model.input, time.taken)
  names(list.all)     <- c('all tracks','most probable track','all possible particles','input parameters','model run time') 
  
  options(warn = oldw)
  
  return(list.all)
}
