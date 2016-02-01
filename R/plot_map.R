#' plot map of model output
#' 
#' plot map of model output
#' @param pr promm model output
#' @export


plot_map <- function(pr){
  
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(pr[[1]],col="white")
  map('world',add=T,col=4)
  for(s in 1:length(unique(pr[[2]]$step))){
    plot(pr[[1]][pr[[1]]$step==unique(pr[[1]]$step)[s],],col=colorRampPalette(c('grey90','grey50'))(nrow(pr[[2]]))[s],
         add=T,pch=19,cex=0.6)
  }
  
  mm2 <-pr[[2]]; lines(mm2$lon,mm2$lat,col=1)
  points(mm2$lon,mm2$lat,cex=1,pch=21,bg=colorRampPalette(c('yellow','darkred'))(nrow(pr[[2]])))
  
  par(mfrow=c(1,1),mar=c(4,4,2,2)) 
  
  
}