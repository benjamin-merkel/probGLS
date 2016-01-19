#' probabalistic movement model for geolocation data
#' 
#' modeling stuff 
#' @param particle.number number of particles for each location cloud used in the model
#' @param bootstrap.number number of iterations
#' @param loess.quartile quartiles for loess.filter (Geolight), if NULL then loess filter is not used
#' @param colony colony longitude and latitude
#' @param tagging.date deployment data
#' @param retrieval.date  retrieval date
#' @param twilight.sd sd around each unsrise or sunset event in min
#' @param range.sun.elev min and max sun angle in degree, resolution of sun elevation angle
#' @param min.lat min latitude where you expect an animal
#' @param max.lat max latitude where you expect an animal
#' @param optimal.speed optimal speed of the animal in m/s
#' @param speed.sd speed sd in m/s
#' @param max.speed.allowed maximum speed allowed in m/s
#' @param sst.sd SST sd in degree C
#' @param max.sst.diff max difference in SST allowed in degree C
#' @param days.around.spring.equinox days before equinox, days after equinox
#' @param days.around.fall.equinox days before equinox, days after equinox
#' @param ice.conc.cutoff min percentage of sea ice in which there are no birds
#' @param plot.it plot each step
#' @param boundary.box min lon, max lon, min lat and max lat of extrem boundary
#' @param med.black.sea if T remove meditearanian and black sea
#' @param baltic.sea if T remove baltic sea
#' @param caspian.sea if T remove caspian sea
#' @param east.west.comp if T apply biotrack east west movement compensation
#' @param slog3 sensor data input
#' @param ho2 sunrise and sunset data input 
#' @param act activity data input
#' @param landmask.location directory location of land mask file
#' @param NOAA.OI.location directory location of NOAA OI V2 files
#' @return A list of modeling results: [1] bootstrapped data, [2] geographic median position, [3] raw data, [4] raw data prepared, [5] all possible particles, [7] daily sensor data, [8] time.taken
#' @export



promm <-  function(particle.number             = 500
                           ,bootstrap.number            = 100
                           ,loess.quartile              = NULL  #quartiles for Geolight, if NULL then loess filter is not used
                           ,colony                      = c(15.15,78.17) # colony longitude and latitude
                           ,tagging.date                = as.POSIXct("2014-06-22")
                           ,retrieval.date              = as.POSIXct("2015-07-05")
                           ,twilight.sd                 = 10           # in min
                           ,range.sun.elev              = c(-7,-1,0.1) # min and max sun angle in degree, resolution of sun elevation angle
                           ,min.lat                     = 40           # min latitude where you expect a bird
                           ,max.lat                     = 90           # max latitude where you expect a bird
                           ,optimal.speed               = 20           # in m/s
                           ,speed.sd                    = 0.2          # in m/s
                           ,max.speed.allowed           = 25           # in m/s
                           ,sst.sd                      = 0.5            # in degree C
                           ,max.sst.diff                = 3            # in degree C
                           ,days.around.spring.equinox  = c(5,5)     # days before equinox, days after equinox
                           ,days.around.fall.equinox    = c(5,5)     # days before equinox, days after equinox
                           ,ice.conc.cutoff             = 0.9          # min percentage of sea ice in which there are no birds
                           ,plot.it                     = F
                           ,boundary.box                = c(-90,70,-55,90) # min lon, max lon, min lat, max lat
                           ,med.black.sea               = T            # if T remove meditearanian and black sea
                           ,baltic.sea                  = T            # if T remove baltic sea
                           ,caspian.sea                 = T            # if T remove caspian sea
                           ,east.west.comp              = T            # if T apply biotrack east west movement compensation
                           ,slog3                       = sensor[sensor$SST.remove==F,]         
                           ,ho2                         = trn 
                           ,act                         = act
                           ,landmask.location           = 'E:/environmental data/SST/lsmask.oisst.v2.nc'
                           ,NOAA.OI.location            = 'E:/environmental data/SST/NOAA OI SST V2'
                           ){

# load libaries ---------
library(raster)
library(rgeos)
library(Imap)
library(sampSurf)
library(maptools)
library(PBSmapping)
library(rgdal)
library(GeoLight)
library(ncdf4)
library(fields)
library(plyr)
library(spatstat)

  
  
start.time <- Sys.time()

# remove all known data -----
ho2   <- ho2  [ho2$tFirst >= as.POSIXct(tagging.date) & ho2$tSecond <= as.POSIXct(retrieval.date),]
slog3 <- slog3[slog3$date >= as.Date(tagging.date)    & slog3$date  <= as.Date(retrieval.date),]
act   <- act  [act$dtime  >= as.POSIXct(tagging.date) & act$dtime   <= as.POSIXct(retrieval.date),]


# define projections-----
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# create spatial colony object-----
col              <- as.data.frame(cbind(colony[1],colony[2]))
colnames(col)    <- c("lon","lat")
coordinates(col) <- cbind(col$lon,col$lat)
proj4string(col) <- proj.latlon
col$dtime        <- tagging.date
col$doy          <- NA
col$jday         <- NA
col$year         <- NA
col$type         <- NA
col$tFirst       <- as.POSIXct(tagging.date)
col$tSecond      <- as.POSIXct(tagging.date)
col$tFirst.err   <- NA
col$tSecond.err  <- NA
col$sun.elev     <- NA
col$step         <- NA
col$bootstrap    <- NA
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
col$wrel         <- NA


# create empty spdf ----
empty.spdf       <- col
empty.spdf$wrel  <- 2


# weighing of speed and sst----
max.w.speed <- dnorm(x=optimal.speed,mean=optimal.speed,sd=speed.sd)

# calc minimum time required to travel average GLS error at maximum allowed speed----
min.time.dry <- 186 / max.speed.allowed/3.6

# remove all twilight outliers----
if(is.null(loess.quartile)==F){
  ho2$loes  <- loessFilter(ho2, plot = T, k = loess.quartile)
  ho2       <- ho2[ho2$loes == T,]
}

# west east movement compensation----
if(east.west.comp==T){
  datetime   <- ho2$dtime
  temp.lon   <- coord(ho2,degElevation=-4,note=F)[,2]
  longitude2 <- -temp.lon
  longitude1 <- -c(NA,temp.lon[1:(length(temp.lon)-1)])
  timedate2  <- as.numeric(datetime)
  timedate1  <- as.numeric(c(NA,datetime[1:(length(datetime)-1)]))
  length     <- abs(as.numeric(difftime(ho2$tSecond,ho2$tFirst,units="sec")))
  m          <- (longitude2-longitude1)/((timedate2-timedate1)) * length / 15 /2
  ewc        <- data.frame(datetime,longitude2,longitude1,timedate2,timedate1,length,m)
  
  ewc$tFirst.corrected  <- ho2$tFirst +(ewc$m)
  ewc$tSecond.corrected <- ho2$tSecond-(ewc$m)
  ewc$tFirst.corrected [is.na(ewc$tFirst.corrected)]  <- ho2$tFirst [is.na(ewc$tFirst.corrected)]
  ewc$tSecond.corrected[is.na(ewc$tSecond.corrected)] <- ho2$tSecond[is.na(ewc$tSecond.corrected)]
  
  ho2$tFirst  <- ewc$tFirst.corrected
  ho2$tSecond <- ewc$tSecond.corrected
}

# subset data to columns of interest-----
ho3       <- data.frame(subset(ho2,select=c(tFirst,tSecond,type,dtime,doy,jday,year,month)))

# duplicate data frame number of sun.elev range-----
sun.elev.steps <- length(seq(range.sun.elev[1],range.sun.elev[2],range.sun.elev[3]))  
ho4            <- data.frame(mapply(rep,ho3,sun.elev.steps))

# duplicate data frame number of particle times----
ho4     <- data.frame(mapply(rep,ho4,ceiling(particle.number/sun.elev.steps)))

# make sure every column is in the right format-----
ho4[,1] <- as.POSIXct(as.numeric(as.character(ho4[,1])),origin="1970-01-01",tz="UTC")
ho4[,2] <- as.POSIXct(as.numeric(as.character(ho4[,2])),origin="1970-01-01",tz="UTC")
ho4[,3] <- as.numeric(as.character(ho4[,3]))
ho4[,4] <- as.POSIXct(as.numeric(as.character(ho4[,4])),origin="1970-01-01",tz="UTC")
ho4[,5] <- as.numeric(as.character(ho4[,5]))
ho4[,6] <- as.numeric(as.character(ho4[,6]))
ho4[,7] <- as.numeric(as.character(ho4[,7]))
ho4[,8] <- as.numeric(as.character(ho4[,8]))
ho4[,9] <- as.numeric(as.character(ho4[,9]))

# add sun elevation angle-----
for(k in seq(range.sun.elev[1],range.sun.elev[2],range.sun.elev[3])){
  rse       <- rep(k,length=(ceiling(particle.number/sun.elev.steps)*length(ho3[,1])))
  if(k == range.sun.elev[1]) rse2<-rse else rse2<-c(rse2,rse)
}
ho4$sun.elev<- rse2

# vary tFirst and tSecond----
ho4$tFirst.er  <- 60*rnorm(length(ho4[,1]), mean = 0, sd = twilight.sd)  # in sec
ho4$tSecond.er <- 60*rnorm(length(ho4[,1]), mean = 0, sd = twilight.sd)  # in sec
ho4$tFirst     <- ho4$tFirst  + ho4$tFirst.er
ho4$tSecond    <- ho4$tSecond + ho4$tSecond.er


# calculate coordinates out of twilight times and sun elevation angle; coord----
new.pos                 <- as.data.frame(coord(ho4,degElevation=ho4$sun.elev,note=F,method='NOAA'))
colnames(new.pos)       <- c('lon','lat')
ho4$lon                 <- new.pos$lon
ho4$lat                 <- new.pos$lat

# define equinox periods----
spring.equinox <- c((80  - days.around.spring.equinox[1]):(80  + days.around.spring.equinox[2]))
fall.equinox   <- c((264 - days.around.fall.equinox[1])  :(264 + days.around.fall.equinox[2]))

# assign random latitudes to equinox periods----
ho4$lat[ho4$doy %in% c(spring.equinox,fall.equinox)] <- sample(seq(min.lat,max.lat,by=0.001),size=length(ho4$lat[ho4$doy %in% c(spring.equinox,fall.equinox)]),replace=T)

# remove sun elevation angle during equinox----
ho4$sun.elev[ho4$doy %in% c(spring.equinox,fall.equinox)] <- NA

# remove all positions outside boundaries-----
ho4    <- ho4[ho4$lon>boundary.box[1] & ho4$lon<boundary.box[2] & ho4$lat>boundary.box[3] & ho4$lat<boundary.box[4],]

# transform to SpatialPointsdataframe and assign unique step to each particle cloud -----
sp6                <- ho4[!is.na(ho4$lat),]
sp6$lon[sp6$lon<0] <- sp6$lon[sp6$lon<0]+360
coordinates(sp6)   <- cbind(sp6$lon,sp6$lat)
proj4string(sp6)   <- CRS(proj.latlon)
sp6$dtime2         <- (sp6$tSecond-sp6$tFirst)/2+sp6$tFirst
sp6                <- sp6[order(sp6$dtime2),]
sp6$loop.step      <- paste(floor(as.numeric(julian(sp6$tFirst))),
                            floor(as.numeric(julian(sp6$tSecond))),
                            sp6$type,sep="-") 

# add unique step id ----
loopstep           <- data.frame(loop.step=unique(sp6$loop.step),step=c(1:length(unique(sp6$loop.step))))
sp6                <- merge(sp6,loopstep,by="loop.step",all=T)

# remove all point clouds smaller then 1/5 of particles ----
jt                 <- data.frame(jt=names(sort(table(sp6$loop.step))),
                                 no=sort(table(sp6$loop.step)))
sp6                <- sp6[sp6$loop.step %in% jt$jt[jt$no>=c(particle.number/5)],]

# remove point clouds with maximum lower then min.lat----
rm.lat           <- data.frame(max.lat=tapply(sp6$lat,sp6$loop.step,max))
rm.lat$loop.step <- rownames(rm.lat)
rm.lat$rm        <- 1
rm.lat$rm[rm.lat$max.lat < min.lat] <- 0
sp7                  <- sp6[sp6$loop.step %in% rm.lat$loop.step[rm.lat$rm==1],]

# remove everything on land-----
landms               <- raster(landmask.location)
sp7$landmask         <- extract(landms,sp7)
sp7                  <- sp7[sp7$landmask==1,]

# remove all point clouds with less than 10 % points outside land ----
jt                   <- data.frame(jt=names(sort(table(sp7$loop.step))),no=sort(table(sp7$loop.step)))
sp8                  <- sp7[sp7$loop.step %in% jt$jt[jt$no>=c(particle.number*0.1)],]

# change longitude from 0/360 to -180/180 ----
grr                  <- as.data.frame(sp8)
grr$lon[grr$lon>180] <- grr$lon[grr$lon>180]-360

# remove baltic sea ---- 
if(baltic.sea==T) grr$lon[grr$lon>14     & grr$lon<33.5 & grr$lat>51.4 & grr$lat<66.2] <- NA

# remove  meditereanian and black sea -----
if(med.black.sea==T) {
  grr$lon[grr$lon>0 & grr$lon<45 & grr$lat>30 & grr$lat<48] <- NA
  grr$lon[grr$lon>(-5) & grr$lon<0.5 & grr$lat>30 & grr$lat<48] <- NA
}

# remove caspian sea ---- 
if(caspian.sea==T) grr$lon[grr$lon>45 & grr$lon<62 & grr$lat>35 & grr$lat<48] <- NA

# create grr object-----
grr              <- grr[!is.na(grr$lon),]
coordinates(grr) <- cbind(grr$lon,grr$lat)
proj4string(grr) <- CRS(proj.latlon)
grr$date         <- grr$dtime

# define date and time as average between tFirst and tSecond-----
grr$dtime        <- as.POSIXct((as.numeric(grr$tSecond)-as.numeric(grr$tFirst))/2,origin=grr$tFirst,tmz='UTC')
grr$jday2        <- floor(grr$jday)
grr              <- grr[order(grr$dtime),]


# create current location list and start with colony and time of tagging-----
col2             <- col
col2$dtime       <- as.POSIXct(tagging.date)
col2$jday        <- as.numeric(julian(col2$dtime))
colt             <- vector("list",length=bootstrap.number)
colt[1:bootstrap.number] <- col2

# loop through each step ----
iter = 0
for(ts in unique(grr$step)){
  step.start <- Sys.time()
  
  if(length(grr$dtime[grr$step == ts])>0){
    
    # select the right time step----
    gr3          <- grr[grr$step==ts,]
    
    # calculate bearing and distance to previous locations----
    gbear        <- lapply (colt,FUN=function(x) bearing(x,gr3))
    gdist        <- lapply (colt,FUN=function(x) spDists(gr3,x,longlat=T)*1000)  
    
    # calculate what fraction of the time the logger was dry  -----        
    fun.time.dry <- function (x) {
      #slo2   <- slog3[slog3$date>=(as.Date(x$dtime)-1) & slog3$date<=(max(as.Date(gr3$dtime))+1),]
      #sumact <- 1-(mean(slo2$wetdry_ratio)-min.time.dry)/24
      #if(sumact > 1) sumact <- 1
      slo2   <- act[act$dtime >= x$tFirst[1] & act$dtime <= max(gr3$tSecond),]
      slo2.time <- abs(as.numeric(difftime(x$tFirst[1],max(gr3$tSecond),units='secs')))
      sumact <- (1 - (sum(slo2$wets0.20)*30 - min.time.dry*60*60)/slo2.time)
      if(sumact > 1) sumact <- 1
      return(sumact)
    }
    gtime.dry    <- lapply (colt,fun.time.dry)  
    
    
    
    # create list with one data frame for each bootstrap iteration-----
    gr2                     <- vector("list",length=bootstrap.number)
    gr2[1:bootstrap.number] <- gr3
    gr2                     <- lapply(gr2,function(x) data.frame(x))
    gr2                     <- mapply(cbind,gr2,gbear=gbear,gdist=gdist,time.dry=gtime.dry, SIMPLIFY = FALSE)
    
    
    # calculate speed and time difference-----
    prev.dtime <- lapply(colt,function(x) x$dtime)
    gr2        <- mapply(cbind,gr2,prev.dtime=prev.dtime, SIMPLIFY = FALSE)
    gspeed     <- lapply (gr2,FUN=function(x) x$gdist/abs(as.numeric(difftime(x$dtime,x$prev.dtime,units="secs"))*x$time.dry))
    time.diff  <- lapply (gr2,FUN=function(x) difftime(x$dtime,x$prev.dtime,units="mins"))
    gr2        <- mapply(cbind,gr2,gspeed=gspeed,time.diff=time.diff, SIMPLIFY = FALSE)
    
    # load needed SST, SST error and ice data-----
    track2       <- data.frame(day=c(as.numeric(as.character(substr(gr3$dtime[1],9,10)))),
                               month=c(as.numeric(as.character(substr(gr3$dtime[1],6,7)))),
                               year=c(gr3$year[1]),
                               lon=c(floor(min(gr3$lon)),ceiling(max(gr3$lon))),
                               lat=c(floor(min(gr3$lat)),ceiling(max(gr3$lat))))
    
    
    fname.sst <- paste(NOAA.OI.location,'/sst.day.mean.' ,year=c(gr3$year[1]),'.v2.nc',sep='')
    fname.err <- paste(NOAA.OI.location,'/sst.day.err.'  ,year=c(gr3$year[1]),'.v2.nc',sep='')
    fname.ice <- paste(NOAA.OI.location,'/icec.day.mean.',year=c(gr3$year[1]),'.v2.nc',sep='')
        
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
    sstdata$ice     <- eii2[,3]
    sstdata$sst.err <- eri2[,3]
    
    # remove sst values in pixels with more than ice.conc.cutoff----
    # except for 2012-8-11 to 2012-8-16 as ice data is fucked up in this period
    if(gr3$date[1]<"2012-08-11" | gr3$date[1]>"2012-08-16") sstdata$V1[sstdata$ice>=ice.conc.cutoff] <-NA
    
    # rasterize sat data-----
    min.lat.sst   <- floor(min(sstdata$Lat))
    max.lat.sst   <- ceiling(max(sstdata$Lat))
    min.lon.sst   <- floor(min(sstdata$Lon))
    max.lon.sst   <- ceiling(max(sstdata$Lon))
    
    if(min.lat.sst<(-90)) min.lat.sst=(-90)
    if(max.lat.sst>( 90)) max.lat.sst=( 90)
    if(min.lon.sst<(-180)) min.lon.sst=(-180)
    if(max.lon.sst>( 180)) max.lon.sst=( 180)
    
    r <- raster(xmn=min.lon.sst,xmx=max.lon.sst,
                ymn=min.lat.sst,ymx=max.lat.sst,
                crs=CRS(proj.latlon),
                resolution=c(0.25,0.25))
    
    coordinates(sstdata)                    <- cbind(sstdata[,2],sstdata[,1])
    proj4string(sstdata)                    <- proj4string(gr3)
    sstdata$V1[is.na(sstdata$V1)]           <- c(-10)
    sstdata$sst.err[is.na(sstdata$sst.err)] <- c(0)
    sstdata$ice[is.na(sstdata$ice)]         <- c(0)
    
    sstd <- rasterize(sstdata,r,'V1')
    errd <- rasterize(sstdata,r,'sst.err')
    iced <- rasterize(sstdata,r,'ice')
    
    
    # extract satellite SSt and tag sst for each particle and calculate difference----
    gr3$sat.ice      <- extract(iced,gr3)
    gr3$sat.sst      <- extract(sstd,gr3)
    gr3$sat.sst.err  <- extract(errd,gr3)
    if(length(slog3$SST[slog3$jday==gr3$jday2[1]])>0)  gr3$tag.sst <- slog3$SST[slog3$jday==gr3$jday2[1]]
    if(length(slog3$SST[slog3$jday==gr3$jday2[1]])==0) gr3$tag.sst <- NA
    gr3$sat.sst[gr3$sat.sst==c(-10)]      <- NA
    gr3$sat.sst.err[is.na(gr3$sat.sst)]   <- NA
    gr3$sst.diff     <- gr3$sat.sst-gr3$tag.sst
    
    
    gr2   <- lapply(gr2,function(x) cbind(x,sat.sst     = gr3$sat.sst,
                                            tag.sst     = gr3$tag.sst,
                                            sst.diff    = gr3$sst.diff,
                                            sat.sst.err = gr3$sat.sst.err,
                                            sat.ice     = gr3$sat.ice))
    gr2   <- lapply(gr2,function(x) x[!is.na(x$sat.sst),])
    
    
    # weigh each particle according to speed and SST
    fun.gsst   <- function(x) (1/(sqrt(2*pi) * (sst.sd+x$sat.sst.err)) * exp(-((x$sst.diff - 0)^2/(2 *(sst.sd+x$sat.sst.err)^2))))/(1/(sqrt(2*pi) * (sst.sd+x$sat.sst.err)) * exp(-((0 - 0)^2/(2 *(sst.sd+x$sat.sst.err)^2))))
    
    fun.gspeed <- function(x) (1/(sqrt(2*pi) * speed.sd) * exp(-((x$gspeed - optimal.speed)^2/(2 *speed.sd^2))))/max.w.speed 
    
    
    wsst   <- lapply(gr2,fun.gsst)
    wspeed <- lapply(gr2,fun.gspeed)
    gr2    <- mapply(cbind,gr2,wspeed=wspeed,wsst=wsst, SIMPLIFY = FALSE)
    gr2    <- lapply(gr2, function(x) {x$wspeed[x$gspeed<=optimal.speed]<-1 
                                       x$wspeed[x$gspeed< 0]<-0
                                       x$wspeed[x$gspeed>max.speed.allowed]<-0 
                                       x$wsst  [x$sst.diff>  max.sst.diff ]<-0 
                                       x$wsst  [x$sst.diff<(-max.sst.diff)]<-0 
                                       return(x)})
    
    
    # if no tag.SST available only use speed weighing----
    if(is.na(gr3$tag.sst[1])==F) gselect  <- lapply(gr2,function(x) x$wspeed*x$wsst)
    if(is.na(gr3$tag.sst[1])==T) gselect  <- lapply(gr2,function(x) x$wspeed)
    gr2      <- mapply(cbind,gr2,grel=gselect, SIMPLIFY = FALSE)
    
    
    # create list for next step particles-----
    new.r        <- vector("list",length=bootstrap.number)
    new.r2       <- vector("list",length=bootstrap.number)
    random.point <- vector("list",length=bootstrap.number)
        
    # choose a specific random point based on weighing of particles
    for(botts in 1:bootstrap.number){
      
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
        new.r[[botts]]$bootstrap    <- botts
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
        new.r[[botts]]   <- empty.spdf
        new.r2[[botts]]  <- subset(data.frame(new.r[[botts]]),select=names(empty.spdf))
        new.r[[botts]]   <- colt[[botts]]
      }
    }
    
    colt<-new.r
    
    # save loop step data in object ----
    iter = iter + 1
    new.r2 <- as.data.frame(Map(f(new.r2), names(new.r2[[1]])))
    #new.r2 <- subset(new.r2,select=names(empty.spdf))
    if(iter==1) newt2 <- new.r2 else newt2 <- rbind(newt2,new.r2)
    
  }

  step.end <- Sys.time()
  step.time <- step.end - step.start
  
  cat(paste(as.Date(gr3$dtime)[1],'  -  ',ts," of ",length(unique(grr$step)),' steps  -  ',round(step.time,2),' sec\n',sep=""))  
}

#backup <- newt2

coordinates(newt2) <- cbind(newt2$lon,newt2$lat)
proj4string(newt2) <- CRS(proj.latlon)
newt2$dtime        <- as.POSIXct(newt2$dtime  ,origin="1970-01-01")
newt2$tFirst       <- as.POSIXct(newt2$tFirst ,origin="1970-01-01")
newt2$tSecond      <- as.POSIXct(newt2$tSecond,origin="1970-01-01")
newt2$jday2        <- as.numeric(julian(newt2$dtime))
newt2              <- newt2[order(newt2$step),]
newt2              <- newt2[newt2$wrel<2,]
newt2$month        <- as.numeric(strftime(newt2$dtime,'%m'))



for(bot in unique(newt2$step)){
  if(bot==min(unique(newt2$step))) {
    newg <- gCentroid(newt2[newt2$step==bot,])
    newg <- SpatialPointsDataFrame(newg,data.frame(newt2[newt2$step==bot,][1,]))
    newg$new.lat <- median(newt2$lat[newt2$step==bot])
    newg$new.lon <- median(newt2$lon[newt2$step==bot])
    newg$new.sst <- median(newt2$sat.sst[newt2$step==bot])
    newg$new.sun <- median(newt2$sun.elev[newt2$step==bot])
  }
  if(bot!=min(unique(newt2$step))) {
    newg2 <-gCentroid(newt2[newt2$step==bot,])
    newg2 <- SpatialPointsDataFrame(newg2,data.frame(newt2[newt2$step==bot,][1,]))
    newg2$new.lat <- median(newt2$lat[newt2$step==bot])
    newg2$new.lon <- median(newt2$lon[newt2$step==bot])
    newg2$new.sst <- median(newt2$sat.sst[newt2$step==bot])
    newg2$new.sun <- median(newt2$sun.elev[newt2$step==bot])
    newg <- spRbind(newg,newg2)
  }
}

# calculate geographic median for each bootstrap cloud----
for(i in unique(newg$step)){
  sf <- data.frame(spDists(newt2[newt2$step==i,1:2],longlat=T),ncol=length(newt2$step[newt2$step==i]))
  sa <- data.frame(sum.dist=rowMeans(sf),bot=seq(1,length(newt2$step[newt2$step==i]),1))
  gmp <- newt2[newt2$step==i,][sa$sum.dist==min(sa$sum.dist),]
  newg$gm.lon[newg$step==i] <- gmp$lon
  newg$gm.lat[newg$step==i] <- gmp$lat
}


newg <- data.frame(newg)
coordinates(newg) <- cbind(newg$gm.lon,newg$gm.lat)
proj4string(newg) <- CRS(proj.latlon)

end.time   <- Sys.time()
time.taken <- end.time - start.time

list.all            <- list(newt2,newg,ho3,ho2,sp6,slog3,time.taken)
names(list.all)     <- c('bootstrapped data (newt2)','median position (newg)',
                         'raw data (ho3)','raw data prepared (ho2)',
                         'all possible particles (sp6)',
                         'daily sensor (slog3)','time.taken') 

return(list.all)
}