#' SST deduction
#' 
#' SST deduction from logger temperature reading
#' @param datetime date time object as POSIXct in UTC
#' @param temp temperature readings 
#' @param temp.range min and max of temperature range
#' @param plot Logical. Plot the results if TRUE
#' @details only works with temp data that are recorded while the logger is submerged in seawater
#' @examples
#'#################################################
#'# example black browed albatross temperature data 
#'#################################################
#' 
#'# sst data ----
#'sen           <- sst_deduction(datetime = BBA_sst$dtime, temp = BBA_sst$temp, temp.range = c(-2,30))
#'
#'summary(sen)
#' @export



sst_deduction <- function(datetime, temp, temp.range = c(-2, 30), plot = T){
  
  #appease R CMD check
  SST <- SST.remove <- NULL
  
  data       <- data.frame(datetime,temp)
  
  if(nrow(data)==0){
    output <- NULL
  } else {
    data$dtime <- as.POSIXct(data$datetime, tz="UTC")
    data$date  <- as.Date(data$datetime)
    
    data$temp[data$temp > temp.range[2]]<-NA
    data$temp[data$temp < temp.range[1]]<-NA
    
    temp           <- data.frame(tapply(data$temp, data$date, median, na.rm=T))
    colnames(temp) <- "SST"
    temp           <- cbind(date=as.Date(rownames(temp)), temp)
    
    temp$sst.step.before  <- temp$SST - c(NA, temp$SST[1:(nrow(temp) - 1)])
    temp$sst.step.after   <- temp$SST - c(temp$SST[2:(nrow(temp))], NA)
    temp$sst.diff         <- temp$sst.step.before * temp$sst.step.after
    temp$sst.diff2        <- temp$sst.diff
    temp$SST.remove       <- F
    
    temp$SST.remove[temp$sst.diff2>12] <- T
    
    output <- subset(temp,select=c(date, SST, SST.remove))
    
    if(plot){
      opar <- par(mfrow=c(1,1),mar=c(3,4,1,1))
      plot(data$dtime,data$temp,pch=19,
           type="o",col=grey(0.8),ylab=expression(temperature ~degree~C),xlab="")
      points(as.POSIXct(output$date),output$SST,pch=as.numeric(as.factor(output$SST.remove))+20,
             type="o",bg=as.numeric(as.factor(output$SST.remove)))
      par(opar)
    }
    
  }
  return(output)
}