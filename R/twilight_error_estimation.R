#' Estimate the assumed error around a twilight event
#' 
#' This function is a visual aid to estimate the assumed error around a twilight event. The determined shape and scale parameters will be used as input for the prob_algorithm function. This function assumes that the error around a twilight event is log normal distributed. 
#' @param shape shape parameter of twilight error (in minutes)
#' @param scale scale parameter of twilight error (in minutes)
#' @param delay estimated delay of observed twilight event compared to true twilight event (in minutes)
#' @details This function can be used to visualise the error distribution around each twilight event which will be used in the prob_algorithm function. The error is assumed to be log normal distributed. In the case of sunrises, only a small error before and a higher possible error with a long tail after the twilight event is assumed (this is reversed in the case of sunsets).
#' @export


twilight_error_estimation <- function(shape = 2.49, scale = 0.94, delay = 0){
  
  dat      <- data.frame(x = seq(0, 1440, 0.01))
  dat$dist <- dlnorm(dat$x, shape, scale)
  dist.max <- max(dat$dist, na.rm=T)
  x0       <- dat$x[dat$dist == dist.max][1] 
  dat$dist <- dat$dist/ dist.max
  x0.5     <- min(dat$x[dat$x > x0 & dat$dist > 0.495 & dat$dist < 0.505], na.rm = T)
  x0.1     <- min(dat$x[dat$x > x0 & dat$dist > 0.095 & dat$dist < 0.105], na.rm = T)
  x0.05    <- min(dat$x[dat$x > x0 & dat$dist > 0.045 & dat$dist < 0.055], na.rm = T)
  
  opar <- par(mfrow = c(1, 1), mar = c(7, 7, 5, 1), cex.lab = 1.5, cex.axis = 1.5, las = 1,  mgp = c(4.8, 2, 1))
  plot (dat$x - x0 + delay, dat$dist, 
        col = 5, lwd = 3, lty = 1, type = "l",ylab = "probability", xlab = "error (in minutes)", xlim = c(-x0 - 5 + delay, x0.05 + 5 + delay),
        yaxs = "i", ylim = c(0, 1.05), axes = F)
  axis(1)
  axis(2)
  abline(v = 0, lty = 2, col = 'grey25', lwd = 4)
  abline(v = c(-x0 + delay), lty = 3, col = "firebrick", lwd = 2)
  lines(c(x0.5 - x0 + delay, x0.5 - x0 + delay),c(0, 0.5), lty = 3, col = "firebrick", lwd = 2)
  lines(c(x0.1 - x0 + delay, x0.1 - x0 + delay),c(0, 0.1), lty = 3, col = "firebrick", lwd = 2)
  mtext(paste("maximum error before twilight event =", round(-x0 + delay, 1),"min",
              "\n50% probability error after twilight event =", round(x0.5 - x0 + delay, 1),"min",
              "\n10% probability error after twilight event =", round(x0.1 - x0 + delay, 1),"min"), font = 3, col = "firebrick", cex = 1.2, line = 0.5)
  par(opar)
  
  output <- c(shape = shape, scale = scale, x0 = -x0 + delay, x0.5 = x0.5 - x0 + delay, x0.1 = x0.1 - x0 + delay)
  return(output)
}