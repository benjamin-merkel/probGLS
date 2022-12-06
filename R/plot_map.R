#' plot map of algorithm output
#' 
#' plot map of algorithm output
#' @param pr probGLS algorithm output
#' @param legend.position position of colour legend. Options are "topleft", "topright", "bottomleft", and "bottomright"
#' @import rnaturalearth
#' @import paletteer
#' @import adehabitatHR
#' @export


plot_map <- function(pr, legend.position = "topleft"){
  
  pr1 <- pr[[1]]
  pr2 <- pr[[2]]
  
  boundary <- as.numeric(unlist(strsplit(as.character(pr[[4]][17,2]),'[ ]')))
  colony   <- as.numeric(unlist(strsplit(as.character(pr[[4]][4,2]),'[ ]')))
  colony   <- data.frame(lon = colony[1], lat = colony[2])
  colony   <- st_as_sf(colony, coords = c('lon','lat'), crs = 4326)
  
  doy      <- as.numeric(strftime(pr2$dtime, "%j"))
  lab      <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  month_colours <- paletteer_d("ggthemes::Classic_Cyclic", n = 12) 
  doy_colours   <- colorRampPalette(month_colours[c(1:12,1)])(366)
  
  
  proj <- paste0("+proj=aeqd +lat_0=",median(st_coordinates(pr1)[,2])," +lon_0=",median(st_coordinates(pr1)[,1])," +units=km")
  sf::sf_use_s2(FALSE)
  
  land <- ne_countries(type = 'countries', scale = 'medium')
  land <- st_transform(st_crop(st_as_sf(land), extent(boundary)), proj)
  
  colony <- st_transform(colony, proj)
  
  pr1 <- st_transform(pr1, proj)
  pr2 <- st_transform(pr2, proj)
  
  pr1 <- as_Spatial(pr1)
  mcp95<- st_as_sf(mcp(pr1[,"step"], percent = 95))
  mcp75<- st_as_sf(mcp(pr1[,"step"], percent = 75))
  mcp50<- st_as_sf(mcp(pr1[,"step"], percent = 50))
  mcp25<- st_as_sf(mcp(pr1[,"step"], percent = 25))
  
  mcp95_all <- st_union(mcp95)
  mcp75_all <- st_union(mcp75)
  mcp50_all <- st_union(mcp50)
  mcp25_all <- st_union(mcp25)
  
  pr_eq <- pr2[is.na(pr2$median.solar.angle),]
  

  # plot ----
  opar <- par(mfrow=c(1,1),mar=c(4,4,0,0))
  plot(mcp95_all, border = grey(0.2), lty=3, lwd= 0.5, col = grey(0.98), ylab="x [km]",xlab="y [km]")
  plot(mcp75_all, border = "transparent", col = grey(0.94), add=T)
  plot(mcp50_all, border = "transparent", col = grey(0.9), add=T)
  plot(mcp25_all, border = "transparent", col = grey(0.86), add=T)
  
  axis(1, las= 1)
  axis(2, las= 1)
  plot(land, add=T, col= grey(0.5), border= grey(0.5), lwd= 0.5)
  
  lines(st_coordinates(pr2), col=1)
  lines(c(st_coordinates(pr2)[1,1], st_coordinates(colony)[,1]), c(st_coordinates(pr2)[1,2], st_coordinates(colony)[,2]),lty=2)
  lines(c(st_coordinates(pr2)[nrow(pr2),1], st_coordinates(colony)[,1]), c(st_coordinates(pr2)[nrow(pr2),2], st_coordinates(colony)[,2]),lty=2)
  
  points(st_coordinates(pr2), cex=pr2$mean.rel_weight * 2, pch = 21, bg = doy_colours[doy])
  points(st_coordinates(pr_eq), cex = pr_eq$mean.rel_weight, pch = 3)
  plot(colony, cex=2.1, pch=23, bg=grey(1), add=T)
  
  
  if(legend.position == "topleft")     lp <- c(0.1,0.3, 0.8,1)
  if(legend.position == "topright")    lp <- c(0.8,1,0.8,1)
  if(legend.position == "bottomleft")  lp <- c(0.1,0.3,0.1,0.3)
  if(legend.position == "bottomright") lp <- c(0.8,1,0.1,0.3)
  
  par(fig = lp, new =T, mar=rep(1,4))
  pie(rep(1,12), labels = lab, col = month_colours, edges = 1000,
      border = "transparent", add = T, clockwise = T)
  
  par(opar) 
}
