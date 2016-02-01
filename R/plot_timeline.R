#' plot timeline of model output
#' 
#' plot timeline of model output
#' @param pr promm model output
#' @param degElevation sun elevation angle for comparison threshold method
#' @export



plot_timeline <- function(pr,degElevation=-3){
  
  
  ho2 <- pr[[2]]
  ho2$tFirst  <- ho2$tFirst  + ho2$tFirst.err
  ho2$tSecond <- ho2$tSecond + ho2$tSecond.err
  ho2$lon <- coord(tFirst=ho2$tFirst,tSecond=ho2$tSecond,type=ho2$type,degElevation = degElevation,note=F)[,1]
  ho2$lat <- coord(tFirst=ho2$tFirst,tSecond=ho2$tSecond,type=ho2$type,degElevation = degElevation,note=F)[,2]
  
  
  poly.frame<-function(data1,data2,prob1,prob2){
    polyf<-data.frame(c(unique(data1[order(data1)]),unique(data1[order(data1,decreasing=T)])),c(tapply(data2,data1,quantile, probs = prob1,na.rm=T),tapply(data2,data1,quantile, probs = prob2,na.rm=T)[order(as.numeric(names(tapply(data2,data1,quantile, probs = prob2,na.rm=T))),decreasing=T)]))
    return(polyf)
  }
  
  par(mfrow=c(3,1),mar=c(0,4,0,0))
  plot(pr[[1]]$jday,pr[[1]]$lat,col='white',lwd=1,type="p",xaxt="n",ylab="Latitude")
  polygon(poly.frame(pr[[1]]$jday,pr[[1]]$lat,0.95,0.05),col=rgb(1,0,0,alpha=0.5) ,border=NA)
  lines(pr[[2]]$jday,pr[[2]]$lat,col='darkred',lwd=1,type="o")
  lines(ho2$jday,ho2$lat,lwd=1,type="o")
  
  plot(pr[[1]]$jday,pr[[1]]$lon,col='white',lwd=1,type="p",xaxt="n",ylab="Longitude")
  polygon(poly.frame(pr[[1]]$jday,pr[[1]]$lon,0.95,0.05),col=rgb(1,0,0,alpha=0.5) ,border=NA)
  lines(pr[[2]]$jday,pr[[2]]$lon,col='darkred',lwd=1,type="o")
  lines(ho2$jday,ho2$lon,lwd=1,type="o")
  
  par(mar=c(2,4,0,0))
  plot  (pr[[1]]$jday,pr[[1]]$sat.sst,type='l',lwd=1,col="white",ylab="SST",xaxt="n")
  polygon(poly.frame(pr[[1]]$jday,pr[[1]]$sat.sst,0.95,0.05),col=rgb(1,0,0,alpha=0.5) ,border=NA)
  points(pr[[2]]$jday,pr[[2]]$median.sat.sst,type='o',lwd=1,col="darkred")
  points(pr[[2]]$jday,pr[[2]]$tag.sst,type='o',lwd=1)
  
  par(mfrow=c(1,1),mar=c(4,4,2,2)) 
  
}