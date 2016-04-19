#' plot map of algorithm output
#' 
#' plot map of algorithm output
#' @param pr probGLS algorithm output
#' @import maps
#' @export


plot_map <- function(pr){
  
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(pr[[1]],col="white")
  for(s in 1:length(unique(pr[[2]]$step))){
    plot(pr[[1]][pr[[1]]$step==unique(pr[[1]]$step)[s],],col=colorRampPalette(c('grey90','grey50'))(nrow(pr[[2]]))[s],
         add=T,pch=19,cex=0.3)
  }
  map('world',add=T,col=4,lwd=0.5)
  
  mm2 <-pr[[2]]
  lines(mm2$lon,mm2$lat,col=1)
  points(mm2$lon,mm2$lat,cex=1,pch=21,bg=colorRampPalette(c('yellow','darkred'))(nrow(pr[[2]])))
  mm3 <- mm2[is.na(mm2$median.sun.elev),]
  points(mm3$lon,mm3$lat,cex=0.7,pch=3)
  
  par(mfrow=c(1,1),mar=c(4,4,2,2)) 
  
  
}