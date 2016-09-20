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
#' @param range.solar min and max of solar angle range in degree
#' @param speed.dry optimal speed, speed standard deviation and max speed allowed if logger is dry in m/s
#' @param speed.wet optimal speed, speed standard deviation and max speed allowed if logger is wet in m/s
#' @param sst.sd SST standard deviation in degree C 
#' @param max.sst.diff max difference in SST allowed in degree C
#' @param days.around.spring.equinox days before the Spring equinox and days after the Spring equinox. The Spring equinox is assumed constant at 20 March.
#' @param days.around.fall.equinox days before the Fall equinox, days after the Fall equinox. The Fall equinox is assumed constant at 22 September.
#' @param ice.conc.cutoff max percentage of sea ice in which the animal is believed to be
#' @param boundary.box min lon, max lon, min lat and max lat of extrem boundary where you expect an animal to be
#' @param med.sea if T classifiy mediterranean sea as land   
#' @param black.sea if T classifiy black sea as land      
#' @param baltic.sea if T classifiy baltic sea as land   
#' @param caspian.sea if T classifiy caspian sea as land   
#' @param land.mask if T animal is only using ocean areas, if F animal is only using land areas, if NULL no land mask used
#' @param east.west.comp if T apply biotrack east west movement compensation (Biotrack manual v11 page 31pp.)
#' @param sensor data.frame with daily SST data deduced from tag temperature readings (sst_deduction ouput), NULLif no SST data is available (SST will not be used)
#' @param trn data.frame containing twilights and at least tFirst, tSecond and type (same as computed by trn_to_dataframe, ipe_to_dataframe or lotek_to_dataframe)
#' @param act data.frame containing wet dry data (e.g. .act file from Biotrack loggers or .deg file from migrate tech loggers), NULL if no wetdry data is available (algorithm will assume that the logger was always dry)
#' @param wetdry.resolution sampling rate of conductivity switch in sec (e.g. MK15 & MK3006 sample every 3 sec)
#' @param backward run algorithm from end to start
#' @param NOAA.OI.location directory location of NOAA OI V2 NCDF files as well as land mask file 'lsmask.oisst.v2.nc' (downloadable from http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html)
#' @return A list with: [1] all positions, [2] geographic median positions, [3] all possible particles, [4] input parameters, [5] model run time. List items 1 to 3 are returned as SpatialPointsDataframe.
#' @details Many weighting parameters can be used. Some others (which are not yet implemented) are: surface air temperature and topography/ bathymetry.
#' @import raster
#' @import ncdf4
#' @import geosphere
#' @import GeoLight
#' @import sp
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
#'trn           <- twilightCalc(BBA_lux$dtime, BBA_lux$lig, ask = FALSE, LightThreshold = 2)
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
#'# http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
#'# and place all into the same folder
#'
#'# Also, download the land mask file: 'lsmask.oisst.v2.nc' from the same directory 
#'# and place it in the same folder as all the other NetCDF files
#'
#'# run algorithm ----
#'pr   <- prob_algorithm(trn                         = trn, 
#'                       sensor                      = sen, 
#'                       act                         = act, 
#'                       tagging.date                = start, 
#'                       retrieval.date              = end, 
#'                       loess.quartile              = NULL, 
#'                       tagging.location            = c(-36.816,-54.316), 
#'                       particle.number             = 1000, 
#'                       iteration.number            = 100,
#'                       sunrise.sd                  = tw,
#'                       sunset.sd                   = tw,
#'                       range.solar                 = c(-7,-1),
#'                       boundary.box                = c(-120,40,-90,0),
#'                       days.around.spring.equinox  = c(0,0), 
#'                       days.around.fall.equinox    = c(0,0),
#'                       speed.dry                   = c(12,6,45),
#'                       speed.wet                   = c(1,1.3,5), 
#'                       sst.sd                      = 0.5,       
#'                       max.sst.diff                = 3,          
#'                       east.west.comp              = TRUE, 
#'                       ice.conc.cutoff             = 1, 
#'                       wetdry.resolution           = 1,
#'                       NOAA.OI.location            = "folder with environmental data and land mask")
#'
#'# plot lat, lon, SST vs time ----
#'plot_timeline(pr,degElevation = NULL)
#'
#'# plot lon vs lat map ----
#'plot_map(pr)
#' @export



prob_algorithm <- function(particle.number      = 2000
                   ,iteration.number            = 60
                   ,loess.quartile              = NULL 
                   ,tagging.location            = c(0,0)
                   ,tagging.date     
                   ,retrieval.date   
                   ,sunrise.sd                  = c(2.49, 0.94, 4.98)         
                   ,sunset.sd                   = c(2.49, 0.94, 4.98)         
                   ,range.solar                 = c(-7,-1)
                   ,speed.wet                   = c(20,0.2,25)
                   ,speed.dry                   = c(20,0.2,25)
                   ,sst.sd                      = 0.5       
                   ,max.sst.diff                = 3          
                   ,days.around.spring.equinox  = c(10,10)   
                   ,days.around.fall.equinox    = c(10,10) 
                   ,ice.conc.cutoff             = 1
                   ,boundary.box                = c(-180,180,-90,90)
                   ,med.sea                     = T        
                   ,black.sea                   = T        
                   ,baltic.sea                  = T      
                   ,caspian.sea                 = T    
                   ,land.mask                   = NULL
                   ,east.west.comp              = T   
                   ,sensor        
                   ,trn     
                   ,act  
                   ,wetdry.resolution           = 30
                   ,backward                    = F
                   ,NOAA.OI.location            = 'E:/environmental data/SST/NOAA OI SST V2'){

start.time <- Sys.time()

#appease R CMD check
tFirst <- tSecond <- type <- dtime <- doy <- jday <- year <- month <- NULL



model.input <- data.frame(parameter=c('particle.number','iteration.number','loess.quartile','tagging.location',
                                     'tagging.date','retrieval.date','sunrise.sd','sunset.sd','range.solar','speed.wet',
                                     'speed.dry','sst.sd','max.sst.diff','days.around.spring.equinox',
                                     'days.around.fall.equinox','ice.conc.cutoff','boundary.box','med.sea','black.sea',
                                     'baltic.sea','caspian.sea','east.west.comp','wetdry.resolution','NOAA.OI.location','backward'),
                          chosen=c(paste(particle.number,collapse=" "),paste(iteration.number,collapse=" "),paste(loess.quartile,collapse=" "),paste(tagging.location,collapse=" "),
                                   paste(tagging.date,collapse=" "),paste(retrieval.date,collapse=" "),paste(sunrise.sd,collapse=" "),paste(sunset.sd,collapse=" "),paste(range.solar,collapse=" "),paste(speed.wet,collapse=" "),
                                   paste(speed.dry,collapse=" "),paste(sst.sd,collapse=" "),paste(max.sst.diff,collapse=" "),paste(days.around.spring.equinox,collapse=" "),
                                   paste(days.around.fall.equinox,collapse=" "),paste(ice.conc.cutoff,collapse=" "),paste(boundary.box,collapse=" "),paste(med.sea,collapse=" "),paste(black.sea,collapse=" "),
                                   paste(baltic.sea,collapse=" "),paste(caspian.sea,collapse=" "),paste(east.west.comp,collapse=" "),paste(wetdry.resolution,collapse=" "),paste(NOAA.OI.location,collapse=" "),
                                   paste(backward,collapse=" ")))

# find land mask file or error ----
landmask.location <- list.files(path=NOAA.OI.location,pattern="lsmask.oisst.v2.nc",recursive=T)
if(length(landmask.location)==0){
  stop(paste('no land mask file found in folder',NOAA.OI.location,sep=' '),call.=F)
}
landmask.location <- paste(NOAA.OI.location,landmask.location,sep='/')[1]



# add date time julian doy etc -----
trn$dtime     <- trn$tFirst+as.numeric(difftime(trn$tSecond,trn$tFirst,units='sec'))/2
trn$doy       <- as.numeric(strftime(trn$dtime, format = "%j"))
trn$month     <- as.numeric(strftime(trn$dtime, format = "%m"))
trn$year      <- as.numeric(strftime(trn$dtime, format = "%Y"))
trn$jday      <- as.numeric(julian(trn$dtime))

# remove outside known data -----
trn    <- trn[trn$tFirst >= as.POSIXct(tagging.date) & trn$tSecond <= as.POSIXct(retrieval.date),]
trn    <- trn[!is.na(trn$tFirst),]
trn    <- trn[!is.na(trn$tSecond),]


if(nrow(trn)==0)  stop('no data points in trn file between selected tagging and retrieval date',call.=F)

# define projections-----
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# create empty global raster for NOAA OISST V2 data set ----
r <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,crs=CRS(proj.latlon),resolution=c(0.25,0.25))

# add function to collapse lists into data.frame ----
f = function(x) function(i) sapply(x, `[[`, i)


# create spatial tagging.location object-----
col              <- as.data.frame(cbind(as.numeric(tagging.location[1]),as.numeric(tagging.location[2])))
colnames(col)    <- c("lon","lat")
coordinates(col) <- cbind(col$lon,col$lat)
proj4string(col) <- proj.latlon
col$dtime        <- as.POSIXct(tagging.date,tz="UTC")
col$doy          <- NA
col$jday         <- NA
col$year         <- NA
col$type         <- NA
col$tFirst       <- as.POSIXct(tagging.date,tz="UTC")
col$tSecond      <- as.POSIXct(tagging.date,tz="UTC")
col$tFirst.err   <- NA
col$tSecond.err  <- NA
col$sun.elev     <- NA
col$step         <- NA
col$iteration    <- NA
col$bearing      <- NA
col$distance     <- NA
col$frac.timedry <- NA
col$speed        <- NA
col$timediff     <- NA
col$sat.ice      <- NA
col$sat.sst      <- NA
col$sat.sst.err  <- NA
col$tag.sst      <- NA
col$sst.diff     <- NA
col$wsst         <- NA
col$wspeed       <- NA
col$wrel         <- 1

# create empty spdf ----
empty.spdf       <- col
empty.spdf$wrel  <- 2

# remove all twilight outliers----
if(!is.null(loess.quartile)){
  trn$loes  <- loessFilter(trn, plot = T, k = loess.quartile)
  trn       <- trn[trn$loes == T,]
}

# west east movement compensation----
if(east.west.comp==T){
  datetime   <- trn$dtime
  temp.lon   <- coord(trn,degElevation= -6,note=F)[,2]
  longitude2 <- -temp.lon
  longitude1 <- -c(NA,temp.lon[1:(length(temp.lon)-1)])
  timedate2  <- as.numeric(datetime)
  timedate1  <- as.numeric(c(NA,datetime[1:(length(datetime)-1)]))
  length     <- abs(as.numeric(difftime(trn$tSecond,trn$tFirst,units="sec")))
  m          <- (longitude2-longitude1)/((timedate2-timedate1)) * length / 15 /2
  ewc        <- data.frame(datetime,longitude2,longitude1,timedate2,timedate1,length,m)
  
  ewc$tFirst.corrected  <- trn$tFirst +(ewc$m)
  ewc$tSecond.corrected <- trn$tSecond-(ewc$m)
  ewc$tFirst.corrected [is.na(ewc$tFirst.corrected)]  <- trn$tFirst [is.na(ewc$tFirst.corrected)]
  ewc$tSecond.corrected[is.na(ewc$tSecond.corrected)] <- trn$tSecond[is.na(ewc$tSecond.corrected)]
  
  trn$tFirst  <- ewc$tFirst.corrected
  trn$tSecond <- ewc$tSecond.corrected
}

# subset data to columns of interest-----
ho3     <- data.frame(subset(trn,select=c(tFirst,tSecond,type,dtime,doy,jday,year,month)))

# duplicate data frame number of particle times----
ho4     <- data.frame(mapply(rep,ho3,particle.number))

# make sure every column is in the right format-----
ho4[,1] <- as.POSIXct(as.numeric(as.character(ho4[,1])),origin="1970-01-01",tz="UTC")
ho4[,2] <- as.POSIXct(as.numeric(as.character(ho4[,2])),origin="1970-01-01",tz="UTC")
ho4[,4] <- as.POSIXct(as.numeric(as.character(ho4[,4])),origin="1970-01-01",tz="UTC")

# assign unique step to each particle cloud -----
ho4$loop.step      <- paste(as.numeric(julian(ho4$tFirst)),as.numeric(julian(ho4$tSecond)),ho4$type,sep="-") 
ho4$step           <- as.numeric(as.factor(ho4$loop.step)) 

# add sun elevation angle-----
sun.elev.steps <- seq(range.solar[1],range.solar[2],0.01)
ho4$sun.elev   <- sample(sun.elev.steps,size=nrow(ho4),replace=T)

# vary tFirst and tSecond----
ho4$tFirst.er [ho4$type==1] <-  60 * (rlnorm(length(ho4[ho4$type==1,1]), meanlog = sunrise.sd[1], sdlog = sunrise.sd[2]) + sunrise.sd[3])
ho4$tSecond.er[ho4$type==1] <- -60 * (rlnorm(length(ho4[ho4$type==1,1]), meanlog = sunset.sd[1],  sdlog = sunset.sd[2])  + sunset.sd[3]) 
ho4$tFirst.er [ho4$type==2] <- -60 * (rlnorm(length(ho4[ho4$type==2,1]), meanlog = sunset.sd[1],  sdlog = sunset.sd[2])  + sunset.sd[3]) 
ho4$tSecond.er[ho4$type==2] <-  60 * (rlnorm(length(ho4[ho4$type==2,1]), meanlog = sunrise.sd[1], sdlog = sunrise.sd[2]) + sunrise.sd[3]) 


ho4$tFirst     <- ho4$tFirst  + ho4$tFirst.er
ho4$tSecond    <- ho4$tSecond + ho4$tSecond.er

# calculate coordinates out of twilight times and sun elevation angle; coord----
new.pos                 <- as.data.frame(coord(ho4,degElevation=ho4$sun.elev,note=F,method='NOAA'))
colnames(new.pos)       <- c('lon','lat')
ho4$lon                 <- new.pos$lon
ho4$lat                 <- new.pos$lat

# define equinox periods----
# doy 79  = 20 March
# doy 265 = 22 September
spring.equinox <- c((79  - days.around.spring.equinox[1]):(79  + days.around.spring.equinox[2]))
fall.equinox   <- c((265 - days.around.fall.equinox[1])  :(265 + days.around.fall.equinox[2]))

# assign random latitudes to equinox periods----
ho4$lat[ho4$doy %in% c(spring.equinox,fall.equinox)] <- sample(seq(boundary.box[3],boundary.box[4],by=0.0001),size=length(ho4$lat[ho4$doy %in% c(spring.equinox,fall.equinox)]),replace=T)

# remove sun elevation angle during equinox----
ho4$sun.elev[ho4$doy %in% c(spring.equinox,fall.equinox)] <- NA

# remove all positions outside boundaries-----
ho4    <- ho4[ho4$lon>boundary.box[1] & ho4$lon<boundary.box[2] & ho4$lat>boundary.box[3] & ho4$lat<boundary.box[4],]

# transform to SpatialPointsdataframe -----
sp6                <- ho4[!is.na(ho4$lat),]

# remove all point clouds smaller then 1/5 of particles ----
jt                 <- data.frame(table(sp6$step))
sp6                <- sp6[sp6$step %in% as.numeric(as.character(jt$Var1[jt$Freq>=c(particle.number/5)])),]

# remove point clouds with maximum lower then min.lat or minimum higher then max.lat----
rm.lat           <- data.frame(max.lat=tapply(sp6$lat,sp6$step,max))
rm.lat$step      <- rownames(rm.lat)

if(nrow(rm.lat)==0){
  stop(paste("No data points inside boundary box. increase boundary box"),call.=F)  
}

rm.lat$rm        <- 1
rm.lat$rm[rm.lat$max.lat < boundary.box[3]] <- 0
sp7              <- sp6[sp6$step %in% rm.lat$step[rm.lat$rm==1],]

rm.lat           <- data.frame(min.lat=tapply(sp7$lat,sp7$step,min))
rm.lat$step      <- rownames(rm.lat)

if(nrow(rm.lat)==0){
  stop(paste("No data points inside boundary box. increase boundary box"),call.=F)  
}

rm.lat$rm        <- 1
rm.lat$rm[rm.lat$min.lat > boundary.box[4]] <- 0
sp7                  <- sp7[sp7$step %in% rm.lat$step[rm.lat$rm==1],]

# create data frame of all particles computed ----
all.particles        <- data.frame(sp7)
#all.particles$lon[all.particles$lon>180] <- all.particles$lon[all.particles$lon>180]-360
coordinates(all.particles) <- cbind(all.particles$lon,all.particles$lat)
proj4string(all.particles) <- CRS(proj.latlon)

# remove everything on land-----
landms               <- rotate(raster(landmask.location))

coordinates(sp7)   <- cbind(sp7$lon,sp7$lat)
proj4string(sp7)   <- CRS(proj.latlon)
sp7$landmask       <- extract(landms,sp7)


# remove baltic sea ---- 
if(baltic.sea==T)  sp7$landmask[sp7$lon>14     & sp7$lon<33.5 & sp7$lat>51.4 & sp7$lat<66.2] <- 0

# remove meditereanian sea -----
if(med.sea==T) {
  sp7$landmask[sp7$lon>=0  & sp7$lon<=27  & sp7$lat>30 & sp7$lat<48] <- 0
  sp7$landmask[sp7$lon>=27 & sp7$lon< 40  & sp7$lat>30 & sp7$lat<40] <- 0
  sp7$landmask[sp7$lon>355 & sp7$lon<=360 & sp7$lat>30 & sp7$lat<42] <- 0
}

# remove black sea ---- 
if(black.sea==T)   sp7$landmask[sp7$lon>27 & sp7$lon<45 & sp7$lat>40 & sp7$lat<48] <- 0

# remove caspian sea ---- 
if(caspian.sea==T) sp7$landmask[sp7$lon>45 & sp7$lon<62 & sp7$lat>35 & sp7$lat<48] <- 0


if(!is.null(land.mask)){
  if(land.mask==T) sp7 <- sp7[sp7$landmask==1,]
  if(land.mask==F) sp7 <- sp7[sp7$landmask==0,]
}

# remove all point clouds with less than 10 % points outside land ----
jt               <- data.frame(table(sp7$step))
grr              <- sp7[sp7$step %in% jt$Var1[jt$Freq>=c(particle.number*0.1)],]

# define date and time as average between tFirst and tSecond-----
grr$dtime        <- as.POSIXct((as.numeric(grr$tSecond)-as.numeric(grr$tFirst))/2,origin=grr$tFirst,tmz='UTC')
grr$jday2        <- floor(grr$jday)
grr$date         <- grr$dtime
grr              <- grr[order(grr$dtime),]

# create current location list and start with tagging.location and time of tagging-----
col2             <- col
if(backward==F) col2$dtime       <- as.POSIXct(tagging.date)
if(backward==T) col2$dtime       <- as.POSIXct(retrieval.date)
col2$jday        <- as.numeric(julian(col2$dtime))
colt             <- vector("list",length=iteration.number)
colt[1:iteration.number] <- col2

# loop through each step ----
iter = 0

if(backward==F) steps <- sort(unique(grr$step))
if(backward==T) steps <- sort(unique(grr$step),decreasing=T)

for(ts in steps){
  step.start <- Sys.time()
  
  if(length(grr$dtime[grr$step == ts])>0){
    
    # select the right time step----
    gr3          <- grr[grr$step==ts,]
    
    # calculate bearing and distance to previous locations----
    gbear        <- lapply (colt,FUN=function(x) bearing(x,gr3))
    gdist        <- lapply (colt,FUN=function(x) spDists(gr3,x,longlat=T)*1000)  
    
    # calculate what fraction of the time the logger was dry  -----        
    fun.time.dry <- function (x) {
      if(!is.null(act)){
      slo2       <- act$wetdry[act$dtime >= min(x$tFirst) & act$dtime <= max(gr3$tSecond)]
      slo2.time  <- abs(as.numeric(difftime(min(x$tFirst), max(gr3$tSecond), units='secs')))
      sumact     <- (1 - sum(slo2) * wetdry.resolution / slo2.time)
      if(sumact > 1) sumact <- 1
      if(sumact < 0) sumact <- 0
      } else {sumact <- 1}
      return(sumact)
    }
    gtime.dry    <- lapply (colt,fun.time.dry)  
    
    # create list with one data frame for each iteration-----
    gr2                     <- vector("list",length=iteration.number)
    gr2[1:iteration.number] <- gr3
    gr2                     <- lapply(gr2,function(x) data.frame(x))
    gr2                     <- mapply(cbind,gr2,gbear=gbear,gdist=gdist,time.dry=gtime.dry, SIMPLIFY = FALSE)
    
    # calculate speed and time difference-----
    prev.dtime <- lapply(colt,function(x) x$dtime)
    gr2        <- mapply(cbind,gr2,prev.dtime=prev.dtime, SIMPLIFY = FALSE)
    gspeed     <- lapply(gr2,function(x) x$gdist / abs(as.numeric(difftime(x$dtime,x$prev.dtime,units="secs"))))
    time.diff  <- lapply(gr2,function(x) difftime(x$dtime,x$prev.dtime,units="mins"))
    gr2        <- mapply(cbind,gr2,gspeed=gspeed,time.diff=time.diff, SIMPLIFY = FALSE)
    
    if(!is.null(sensor)){
      
      sensor$jday  <- as.numeric(julian(sensor$date))
      
      
      # load needed SST, SST error and ice data-----
      track2       <- data.frame(day   = c(as.numeric(as.character(substr(gr3$dtime[1],9,10)))),
                                 month = c(as.numeric(as.character(substr(gr3$dtime[1],6,7)))),
                                 year  = c(gr3$year[1]),
                                 lon   = c(floor(min(gr3$lon)),ceiling(max(gr3$lon))),
                                 lat   = c(floor(min(gr3$lat)),ceiling(max(gr3$lat))))
      
      
      fname.sst <- paste(NOAA.OI.location,'/sst.day.mean.' ,year=track2$year[1],'.v2.nc',sep='')
      fname.err <- paste(NOAA.OI.location,'/sst.day.err.'  ,year=track2$year[1],'.v2.nc',sep='')
      fname.ice <- paste(NOAA.OI.location,'/icec.day.mean.',year=track2$year[1],'.v2.nc',sep='')
      
      if((min(track2$lon)*max(track2$lon))>=0){
        if(min(track2$lon)<0) track2$lon[track2$lon==0]<-360
        track2$lon[track2$lon<0] <- 360+track2$lon[track2$lon<0]
        
        
        eoi <- load.NOAA.OISST.V2 (fname = fname.sst,
                                  lsmask= landmask.location,
                                  lonW  = min(track2$lon),
                                  lonE  = max(track2$lon),
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='sst')
        
        eii <- load.NOAA.OISST.V2 (fname = fname.ice,
                                  lsmask= landmask.location,
                                  lonW  = min(track2$lon),
                                  lonE  = max(track2$lon),
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='icec')
        
        eri <- load.NOAA.OISST.V2 (fname = fname.err,
                                  lsmask= landmask.location,
                                  lonW  = min(track2$lon),
                                  lonE  = max(track2$lon),
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='err')
        
        
        
        eoi2  <-  as.data.frame(as.table(as.matrix(eoi[,,1])))
        colnames(eoi2) <- c("Lat","Long",'V1') 
        eoi2$Lat  <- as.numeric(as.character(eoi2$Lat))
        eoi2$Long <- as.numeric(as.character(eoi2$Long))
        eoi2$Long[eoi2$Long>180] <- eoi2$Long[eoi2$Long>180]-360
        
        eii2  <-  as.data.frame(as.table(as.matrix(eii[,,1])))
        colnames(eii2) <- c("Lat","Long",'V1') 
        eii2$Lat  <- as.numeric(as.character(eii2$Lat))
        eii2$Long <- as.numeric(as.character(eii2$Long))
        eii2$Long[eii2$Long>180] <- eii2$Long[eii2$Long>180]-360
        
        eri2  <-  as.data.frame(as.table(as.matrix(eri[,,1])))
        colnames(eri2) <- c("Lat","Long",'V1') 
        eri2$Lat  <- as.numeric(as.character(eri2$Lat))
        eri2$Long <- as.numeric(as.character(eri2$Long))
        eri2$Long[eri2$Long>180] <- eri2$Long[eri2$Long>180]-360
        
        
      } else {
        track2$lon[track2$lon<0] <- 360+track2$lon[track2$lon<0]
        eoi <- load.NOAA.OISST.V2 (fname = fname.sst,
                                  lsmask= landmask.location,
                                  lonW  = 0,
                                  lonE  = min(track2$lon),
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")))
        
        eii <- load.NOAA.OISST.V2 (fname = fname.ice,
                                  lsmask= landmask.location,
                                  lonW  = 0,
                                  lonE  = min(track2$lon),
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='icec')
        
        eri <- load.NOAA.OISST.V2 (fname = fname.err,
                                  lsmask= landmask.location,
                                  lonW  = 0,
                                  lonE  = min(track2$lon),
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='err')
        
        
        eoi2  <-  as.data.frame(as.table(as.matrix(eoi[,,1])))
        colnames(eoi2) <- c("Lat","Long",'V1') 
        eii2  <-  as.data.frame(as.table(as.matrix(eii[,,1])))
        colnames(eii2) <- c("Lat","Long",'V1') 
        eri2  <-  as.data.frame(as.table(as.matrix(eri[,,1])))
        colnames(eri2) <- c("Lat","Long",'V1') 
        
        eoi <- load.NOAA.OISST.V2 (fname = fname.sst,
                                  lsmask= landmask.location,
                                  lonW  = max(track2$lon),
                                  lonE  = 360,
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")))
        
        eii <- load.NOAA.OISST.V2 (fname = fname.ice,
                                  lsmask= landmask.location,
                                  lonW  = max(track2$lon),
                                  lonE  = 360,
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='icec')
        
        eri <- load.NOAA.OISST.V2 (fname = fname.err,
                                  lsmask= landmask.location,
                                  lonW  = max(track2$lon),
                                  lonE  = 360,
                                  latS  = min(track2$lat),
                                  latN  = max(track2$lat),
                                  date1 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  date2 = as.Date(paste(track2$year[1],track2$month[1],track2$day[1],sep="-")),
                                  extract.value='err')
        
        
        eoi3  <-  as.data.frame(as.table(as.matrix(eoi[,,1])))
        colnames(eoi3) <- c("Lat","Long",'V1') 
        eii3  <-  as.data.frame(as.table(as.matrix(eii[,,1])))
        colnames(eii3) <- c("Lat","Long",'V1') 
        eri3  <-  as.data.frame(as.table(as.matrix(eri[,,1])))
        colnames(eri3) <- c("Lat","Long",'V1') 
        
        eoi2<-rbind(eoi2,eoi3)
        eoi2$Lat  <- as.numeric(as.character(eoi2$Lat))
        eoi2$Long <- as.numeric(as.character(eoi2$Long))
        eoi2$Long[eoi2$Long>180] <- eoi2$Long[eoi2$Long>180]-360
        
        eii2<-rbind(eii2,eii3)
        eii2$Lat  <- as.numeric(as.character(eii2$Lat))
        eii2$Long <- as.numeric(as.character(eii2$Long))
        eii2$Long[eii2$Long>180] <- eii2$Long[eii2$Long>180]-360
        
        eri2<-rbind(eri2,eri3)
        eri2$Lat  <- as.numeric(as.character(eri2$Lat))
        eri2$Long <- as.numeric(as.character(eri2$Long))
        eri2$Long[eri2$Long>180] <- eri2$Long[eri2$Long>180]-360
      }
      
      sstdata         <- eoi2
      sstdata$sst     <- sstdata$V1
      sstdata$ice     <- eii2[,3]
      sstdata$sst.err <- eri2[,3]
      
      # remove sst values in pixels with more than ice.conc.cutoff----
      # except for 2012-8-11 to 2012-8-16 as ice data is fucked up in this period
      if(as.Date(gr3$date[1])<as.Date("2012-08-11") | as.Date(gr3$date[1])>as.Date("2012-08-16")) sstdata$sst[sstdata$ice > ice.conc.cutoff] <-NA
      
      # rasterize sat data-----
      coordinates(sstdata)                    <- cbind(sstdata[,2],sstdata[,1])
      proj4string(sstdata)                    <- CRS(proj.latlon)
      sstdata$sst[is.na(sstdata$sst)]           <- c(-10)
      sstdata$sst.err[is.na(sstdata$sst.err)] <- c(0)
      sstdata$ice[is.na(sstdata$ice)]         <- c(0)
      
      sstd <- rasterize(sstdata,r,'sst')
      errd <- rasterize(sstdata,r,'sst.err')
      iced <- rasterize(sstdata,r,'ice')
      
      
      # extract satellite SSt and tag sst for each particle and calculate difference----
      gr3$sat.ice      <- extract(iced,gr3)
      gr3$sat.sst      <- extract(sstd,gr3)
      gr3$sat.sst.err  <- extract(errd,gr3)
      if(length(sensor$SST[sensor$jday==gr3$jday2[1]])>0)  gr3$tag.sst <- sensor$SST[sensor$jday==gr3$jday2[1]]
      if(length(sensor$SST[sensor$jday==gr3$jday2[1]])==0) gr3$tag.sst <- NA
      gr3$sat.sst[gr3$sat.sst==c(-10)]      <- NA
      gr3$sat.sst.err[is.na(gr3$sat.sst)]   <- NA
      
    } 
    if( is.null(sensor)){
      gr3$sat.ice      <- 0
      gr3$sat.sst      <- 0
      gr3$sat.sst.err  <- 0
      gr3$tag.sst      <- NA
    }
    
    gr3$sst.diff     <- gr3$sat.sst-gr3$tag.sst
    
    
    gr2   <- lapply(gr2,function(x) cbind(x,sat.sst     = gr3$sat.sst,
                                            tag.sst     = gr3$tag.sst,
                                            sst.diff    = gr3$sst.diff,
                                            sat.sst.err = gr3$sat.sst.err,
                                            sat.ice     = gr3$sat.ice))
    if(!is.null(land.mask)) if(land.mask==T) gr2   <- lapply(gr2,function(x) x[!is.na(x$sat.sst),])
    
    
    # weigh each particle according to speed and SST
    
    fun.gspeed <- function(x) {      dnorm(x$gspeed,
                                           mean = c((speed.dry[1]*x$time.dry)+(speed.wet[1]*(1-x$time.dry))),
                                           sd   = c((speed.dry[2]*x$time.dry)+(speed.wet[2]*(1-x$time.dry))))/
                                 max(dnorm(c((speed.dry[1]*x$time.dry)+(speed.wet[1]*(1-x$time.dry))), 
                                           mean = c((speed.dry[1]*x$time.dry)+(speed.wet[1]*(1-x$time.dry))),
                                           sd   = c((speed.dry[2]*x$time.dry)+(speed.wet[2]*(1-x$time.dry)))),na.rm=T)}
    
    wspeed <- lapply(gr2,fun.gspeed)
    
    
    fun.gsst   <- function(x) {   dnorm(x$sst.diff, mean = 0, sd = sst.sd+x$sat.sst.err)/
                              max(dnorm(0     , mean = 0, sd = sst.sd+x$sat.sst.err),na.rm=T)}
    
    wsst <- lapply(gr2,fun.gsst)
    
    
    gr2    <- mapply(cbind, gr2, wspeed = wspeed, wsst = wsst, SIMPLIFY = FALSE)
    
    
    gr2    <- lapply(gr2, function(x) {x$wspeed[x$gspeed<= c(speed.dry[1]*x$time.dry+speed.wet[1]*(1-x$time.dry))]<-1 
                                       x$wspeed[x$gspeed< 0]<-0
                                       x$wspeed[x$gspeed> c(speed.dry[3]*x$time.dry+speed.wet[3]*(1-x$time.dry))]<-0 
                                       x$wspeed[x$sat.ice >ice.conc.cutoff]<- 0 
                                       x$wsst  [x$sst.diff>  max.sst.diff ]<- 0 
                                       x$wsst  [x$sst.diff<(-max.sst.diff)]<- 0 
                                       x$wsst  [is.na(x$wsst)]<- 0 
                                       return(x)})
    
    
    # if no tag.SST available only use speed weighing----
    if(!is.na(gr3$tag.sst[1])) gselect  <- lapply(gr2,function(x) x$wspeed * x$wsst)
    if( is.na(gr3$tag.sst[1])) gselect  <- lapply(gr2,function(x) x$wspeed)
    gr2      <- mapply(cbind,gr2,grel=gselect, SIMPLIFY = FALSE)
    
    
    # create list for next step particles-----
    new.r        <- vector("list",length=iteration.number)
    new.r2       <- vector("list",length=iteration.number)
    random.point <- vector("list",length=iteration.number)
        
    # choose a specific random point based on weighing of particles
    for(botts in 1:iteration.number){
      
      # if all weights are 0 and/or NA jump this position
      if(length(gr2[[botts]]$grel[gr2[[botts]]$grel>0 & !is.na(gr2[[botts]]$grel)])>0){
        
        random.point[[botts]]       <- sample(nrow(gr2[[botts]]),size=1,prob=gr2[[botts]]$grel)
        new.r[[botts]]              <- data.frame(destPoint(colt[[botts]],gr2[[botts]]$gbear[random.point[[botts]]],gr2[[botts]]$gdist[random.point[[botts]]]))
        
        coordinates(new.r[[botts]]) <- new.r[[botts]]
        proj4string(new.r[[botts]]) <- CRS(proj.latlon)
        new.r[[botts]]$dtime        <- gr2[[botts]]$dtime[random.point[[botts]]]
        new.r[[botts]]$doy          <- gr2[[botts]]$doy[random.point[[botts]]]
        new.r[[botts]]$jday         <- gr2[[botts]]$jday[random.point[[botts]]]
        new.r[[botts]]$year         <- gr2[[botts]]$year[random.point[[botts]]]
        new.r[[botts]]$type         <- gr2[[botts]]$type[random.point[[botts]]]
        new.r[[botts]]$tFirst       <- gr2[[botts]]$tFirst[random.point[[botts]]]
        new.r[[botts]]$tSecond      <- gr2[[botts]]$tSecond[random.point[[botts]]]
        new.r[[botts]]$tFirst.err   <- gr2[[botts]]$tFirst.er[random.point[[botts]]]
        new.r[[botts]]$tSecond.err  <- gr2[[botts]]$tSecond.er[random.point[[botts]]]
        new.r[[botts]]$sun.elev     <- gr2[[botts]]$sun.elev[random.point[[botts]]]
        new.r[[botts]]$sex          <- gr2[[botts]]$sex[random.point[[botts]]]
        new.r[[botts]]$morph        <- gr2[[botts]]$morp[random.point[[botts]]]
        new.r[[botts]]$step         <- ts
        new.r[[botts]]$iteration    <- botts
        new.r[[botts]]$bearing      <- gr2[[botts]]$gbear[random.point[[botts]]]
        new.r[[botts]]$distance     <- gr2[[botts]]$gdist[random.point[[botts]]]
        new.r[[botts]]$frac.timedry <- gr2[[botts]]$time.dry[random.point[[botts]]]
        new.r[[botts]]$speed        <- gr2[[botts]]$gspeed[random.point[[botts]]]
        new.r[[botts]]$timediff     <- gr2[[botts]]$time.diff[random.point[[botts]]]
        new.r[[botts]]$sat.ice      <- gr2[[botts]]$sat.ice[random.point[[botts]]]
        new.r[[botts]]$sat.sst      <- gr2[[botts]]$sat.sst[random.point[[botts]]]
        new.r[[botts]]$sat.sst.err  <- gr2[[botts]]$sat.sst.err[random.point[[botts]]]
        new.r[[botts]]$tag.sst      <- gr2[[botts]]$tag.sst[random.point[[botts]]]
        new.r[[botts]]$sst.diff     <- gr2[[botts]]$sst.diff[random.point[[botts]]]
        new.r[[botts]]$wsst         <- gr2[[botts]]$wsst[random.point[[botts]]]
        new.r[[botts]]$wspeed       <- gr2[[botts]]$wspeed[random.point[[botts]]]
        new.r[[botts]]$wrel         <- gr2[[botts]]$grel[random.point[[botts]]]
      
        new.r2[[botts]]  <- subset(data.frame(new.r[[botts]]),select=names(empty.spdf))
      } else {
        new.r [[botts]]  <- empty.spdf
        new.r2[[botts]]  <- subset(data.frame(new.r[[botts]]),select=names(empty.spdf))
        new.r [[botts]]  <- colt[[botts]]
      }
    }
    
    colt<-new.r
    
    # save loop step data in object ----
    iter = iter + 1
    new.r2 <- as.data.frame(Map(f(new.r2), names(new.r2[[1]])))
    
    if(iter==1) newt2 <- new.r2 else newt2 <- rbind(newt2,new.r2)
    
  }

  step.end  <- Sys.time()
  step.time <- step.end - step.start
  
  cat(paste(as.Date(gr3$dtime)[1],'  -  ',iter," of ",length(unique(grr$step)),' steps  -  ',round(step.time,2),' sec',sep=""),'\r')  
}


# remove all empty steps
newt2              <- newt2[newt2$wrel<=1,]
newt2              <- as.data.frame(newt2)


#remove NAs in lat and lon
newt2$lon <- as.numeric(newt2$lon)
newt2$lat <- as.numeric(newt2$lat)
newt2              <- newt2[!is.na(newt2$lon),]
newt2              <- newt2[!is.na(newt2$lat),]

if(is.null(sensor)) {
  newt2$sat.ice     <- NA
  newt2$sat.sst     <- NA
  newt2$sat.sst.err <- NA
  newt2$sst.diff    <- NA
}

#coordinates(newt2) <- c(newt2$lon,newt2$lat)
coordinates(newt2) <- c('lon','lat')
proj4string(newt2) <- CRS(proj.latlon)
newt2$dtime        <- as.POSIXct(newt2$dtime  ,origin="1970-01-01",tz="UTC")
newt2$tFirst       <- as.POSIXct(newt2$tFirst ,origin="1970-01-01",tz="UTC")
newt2$tSecond      <- as.POSIXct(newt2$tSecond,origin="1970-01-01",tz="UTC")
newt2$jday2        <- as.numeric(julian(newt2$dtime))
newt2              <- newt2[order(newt2$step),]
newt2$month        <- as.numeric(strftime(newt2$dtime,'%m'))

# calculate geographic median for each particle cloud----
for(i in unique(newt2$step)){
  sf                  <- data.frame(spDists(newt2[newt2$step==i,1:2],longlat=T),ncol=length(newt2$step[newt2$step==i]))
  sa                  <- data.frame(sum.dist=rowMeans(sf),bot=seq(1,length(newt2$step[newt2$step==i]),1))
  gmp                 <- newt2[newt2$step==i,] [sa$sum.dist==min(sa$sum.dist),] [1,]
  gmp$median.sat.sst  <- median(newt2$sat.sst[newt2$step==i])
  gmp$median.sun.elev <- median(newt2$sun.elev[newt2$step==i])
  gmp$median.wrel     <- median(newt2$wrel[newt2$step==i])
  
  if(i == unique(newt2$step)[1]) newg <- gmp else newg <- rbind(newg,gmp)
}

end.time   <- Sys.time()
time.taken <- abs(difftime(end.time,start.time,units="mins"))

cat(paste('algorithm run time:',round(as.numeric(time.taken),1),'min',sep=" "),'\n')

list.all            <- list(newt2,newg,all.particles,model.input,time.taken)
names(list.all)     <- c('all tracks','most probable track','all possible particles','input parameters','model run time') 

return(list.all)
}