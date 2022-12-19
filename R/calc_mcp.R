#' Calculate MCPs
#' 
#' This function calculates minimum convex polygons (MCP). It is based on the mcp function from adehabitatHR.
#' It does not calculate area and does not require a minimum number of points. It is independent of adehabitatHR and sp, but relies on sf.
#' @param xy SF POINTS data.frame with one column denoting polygon ID.
#' @param percent Percentage considered to determine MCP
#' @param prob1 lower quantile considered
#' @param prob2 upper quantile considered
#' @export

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
    while(nrow(xy.bord) < 4) xy.bord  <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
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

