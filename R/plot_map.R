#' plot map of algorithm output
#' 
#' plot map of algorithm output
#' @param pr probGLS algorithm output
#' @param legend.position position of colour legend. Options are "topleft", "topright", "bottomleft", and "bottomright"
#' @import rnaturalearth
#' @import paletteer
#' @details This function plots the most probable path constructed using prob_algorithm and associated uncertainty. 
#' @details The locations are colour coded by season. First and last location are connected by stippled lines to the capture location (assuming that it is also the recapture location). Dotted shapes denote 95% minimum convex polygon (MCP) of all estimated particles in the track. Increasing grey-scaled shapes illustrate 95, 75, 50 and 25% MCPs of all particles.
#' @importFrom paletteer paletteer_d
#' @importFrom rnaturalearth ne_countries
#' @export

plot_map <- function(pr, legend.position = "topleft"){
  
  x1 <- pr[[1]]
  x2 <- pr[[2]]
  
  boundary <- as.numeric(unlist(strsplit(as.character(pr[[4]][pr[[4]]$parameter=="boundary.box",2]),'[ ]')))
  colony   <- as.numeric(unlist(strsplit(as.character(pr[[4]][pr[[4]]$parameter=="tagging.location",2]),'[ ]')))
  colony   <- data.frame(lon = colony[1], lat = colony[2])
  colony   <- st_as_sf(colony, coords = c('lon','lat'), crs = 4326)
  
  doy      <- as.numeric(strftime(x2$dtime, "%j"))
  lab      <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  month_colours <- paletteer_d("ggthemes::Classic_Cyclic", n = 12) 
  doy_colours   <- colorRampPalette(month_colours[c(1:12,1)])(366)
  
  
  proj <- paste0("+proj=aeqd +lat_0=",median(st_coordinates(x1)[,2])," +lon_0=",median(st_coordinates(x1)[,1])," +units=km")
  sf::sf_use_s2(FALSE)
  
  land <- ne_countries(type = 'countries', scale = 'medium')
  land <- st_transform(st_crop(st_as_sf(land), 
                               xmin = boundary[1], 
                               xmax = boundary[2], 
                               ymin = boundary[3], 
                               ymax = boundary[4]), proj)
  
  colony <- st_transform(colony, proj)
  
  x1 <- st_transform(x1, proj)
  x2 <- st_transform(x2, proj)
  
  mcp95<- calc_mcp(x1[,"step"], percent = 95)
  mcp75<- calc_mcp(x1[,"step"], percent = 75)
  mcp50<- calc_mcp(x1[,"step"], percent = 50)
  mcp25<- calc_mcp(x1[,"step"], percent = 25)
  
  mcp95_all <- st_union(mcp95)
  mcp75_all <- st_union(mcp75)
  mcp50_all <- st_union(mcp50)
  mcp25_all <- st_union(mcp25)
  
  pr_eq <- x2[is.na(x2$median.solar.angle),]
  

  # plot ----
  opar <- par(mfrow=c(1,1),mar=c(4,4,0,0))
  plot(mcp95_all, border = grey(0.2), lty=3, lwd= 0.5, col = grey(0.98), ylab="x [km]",xlab="y [km]")
  plot(mcp75_all, border = "transparent", col = grey(0.94), add=T)
  plot(mcp50_all, border = "transparent", col = grey(0.9), add=T)
  plot(mcp25_all, border = "transparent", col = grey(0.86), add=T)
  
  axis(1, las= 1)
  axis(2, las= 1)
  plot(land, add=T, col= grey(0.5), border= grey(0.5), lwd= 0.5)
  
  lines(st_coordinates(x2), col=1)
  lines(c(st_coordinates(x2)[1,1], st_coordinates(colony)[,1]), c(st_coordinates(x2)[1,2], st_coordinates(colony)[,2]),lty=2)
  lines(c(st_coordinates(x2)[nrow(x2),1], st_coordinates(colony)[,1]), c(st_coordinates(x2)[nrow(x2),2], st_coordinates(colony)[,2]),lty=2)
  
  points(st_coordinates(x2), cex=x2$mean.rel_weight * 2, pch = 21, bg = doy_colours[doy])
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
