#' plot timeline of algorithm output
#' 
#' plot timeline of algorithm output
#' @param pr probGLS algorithm output
#' @param solar.angle sun elevation angle for comparison threshold method, if NULL no threshold positions are plotted
#' @param center.longitude around which longitude should the plot be centered? Only 0 and 180 implemented
#' @export



plot_timeline <- function(pr, solar.angle = NULL, center.longitude = 0){
  
  if(as.character(pr[[4]]$chosen[pr[[4]]$parameter=="sensor.data"])=="TRUE") SST=T else SST=F
  
  if(center.longitude==180){
    x1 <- data.frame(pr[[1]])
    x1$lon[x1$lon<0] <- x1$lon[x1$lon<0]+360
    x2 <- data.frame(pr[[2]])
    x2$lon[x2$lon<0] <- x2$lon[x2$lon<0]+360
    if(!is.null(solar.angle)) {
      x3 <-trn
      x3$lons[x3$lons<0] <- x3$lons[x3$lons<0]+360
    }
    long.label = "Longitude [0 to 360]"
  } else {
    x1 <- pr[[1]]
    x2 <- pr[[2]]
    if(!is.null(solar.angle)) x3 <-trn
    long.label = "Longitude [-180 to 180]"
  }
  x1$lat <- st_coordinates(x1)[,2]
  x1$lon <- st_coordinates(x1)[,1]
  x2$lat <- st_coordinates(x2)[,2]
  x2$lon <- st_coordinates(x2)[,1]
  x2$jday<- as.numeric(julian(x2$dtime))
  
  if(!is.null(solar.angle)){
    trn <- x2
    trn$jday<- as.numeric(julian(trn$dtime))
    pos <- thresholdEstimate(trise  = ifelse(trn$type == 1, trn$tFirst, trn$tSecond),
                             tset   = ifelse(trn$type == 1, trn$tSecond, trn$tFirst),
                             zenith = 90 - solar.angle,
                             tol    = 0.08)
    
    trn$lons <- pos[,1]
    trn$lats <- pos[,2]
  }
  
  
  if(SST==T){
    
    opar <- par(mfrow=c(3,1),mar=c(0,4,0,0),oma=c(2,0,0,0))
    
    plot(x1$jday,x1$lat,col='white',xaxt="n",ylab="Latitude")
    abline(v=x1$jday[is.na(x1$solar.angle)],lty=3,col=grey(0.9),lwd = 0.5)
    polygon(poly_frame(x1$jday,x1$lat,0.75,0.25),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    polygon(poly_frame(x1$jday,x1$lat,0.95,0.05),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    lines(x2$jday,x2$lat,col='darkred',lwd=1,type="o",cex=1)
    if(!is.null(solar.angle)) lines(trn$jday,trn$lats,lwd=1,type="o",cex=1)
    
    plot(x1$jday,x1$lon,col='white',xaxt="n",ylab=long.label)
    polygon(poly_frame(x1$jday,x1$lon,0.75,0.25),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    polygon(poly_frame(x1$jday,x1$lon,0.95,0.05),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    lines(x2$jday,x2$lon,col='darkred',lwd=1,type="o",cex=1)
    if(!is.null(solar.angle)) lines(trn$jday,trn$lons,lwd=1,type="o",cex=1)
    
    
    plot  (x1$jday,x1$sat.sst,col="white",ylab="SST",xaxt="n")
    polygon(poly_frame(x1$jday,x1$sat.sst,0.75,0.25),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    polygon(poly_frame(x1$jday,x1$sat.sst,0.95,0.05),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    points(x2$jday,x2$median.sat.sst,type='o',lwd=1,col="darkred",cex=1)
    points(x2$jday,x2$tag.sst,type='o',lwd=1,cex=1)
    axis(1,at=floor(x2$jday),labels=as.Date(floor(x2$jday),origin="1970-01-01"))
  }
  if(SST==F){
    opar <- par(mfrow=c(2,1),mar=c(0,4,0,0),oma=c(2,0,0,0))
    
    plot(x1$jday,x1$lat,col='white',xaxt="n",ylab="Latitude")
    abline(v=x1$jday[is.na(x1$solar.angle)],lty=3,col=grey(0.9),lwd = 0.5)
    polygon(poly_frame(x1$jday,x1$lat,0.75,0.25),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    polygon(poly_frame(x1$jday,x1$lat,0.95,0.05),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    lines(x2$jday,x2$lat,col='darkred',lwd=1,type="o",cex=1)
    if(!is.null(solar.angle)) lines(trn$jday,trn$lats,lwd=1,type="o",cex=1)
    
    par(mar=c(2,4,0,0))
    plot(x1$jday,x1$lon,col='white',xaxt="n",ylab=long.label)
    polygon(poly_frame(x1$jday,x1$lon,0.75,0.25),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    polygon(poly_frame(x1$jday,x1$lon,0.95,0.05),col=rgb(1,0,0,alpha=0.3) ,border=NA)
    lines(x2$jday,x2$lon,col='darkred',lwd=1,type="o",cex=1)
    if(!is.null(solar.angle)) lines(x3$jday,x3$lons,lwd=1,type="o",cex=1)
    axis(1,at=floor(x2$jday),labels=as.Date(floor(x2$jday),origin="1970-01-01"))
  }
  par(opar) 
  
}
