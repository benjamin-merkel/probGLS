#' SST deduction
#' 
#' SST deduction from logger temperature reading
#' @param datetime date time object as POSIXct in UTC
#' @param temp temperature readings 
#' @param temp.range min and max of temperature range
#' @export



sst_deduction <- function(datetime,temp,temp.range=c(-2,30)){
  
  data       <- data.frame(datetime,temp)
  data$dtime <- as.POSIXct(data$datetime,tz="UTC")
  data$date  <- as.Date(data$datetime)
  
  data$temp[data$temp > temp.range[2]]<-NA
  data$temp[data$temp < temp.range[1]]<-NA
  
  temp           <- data.frame(tapply(data$temp,data$date,median,na.rm=T))
  colnames(temp) <- "SST"
  temp           <- cbind(date=as.Date(rownames(temp)),temp)
  
  temp$doy              <- as.numeric(strftime(temp$date, format = "%j"))
  temp$month            <- as.numeric(strftime(temp$date, format = "%m"))
  temp$year             <- as.numeric(strftime(temp$date, format = "%Y"))
  temp$jday             <- as.numeric(julian(temp$date))
  temp$sst.step.before  <- temp$SST-c(NA,temp$SST[1:(nrow(temp)-1)])
  temp$sst.step.after   <- temp$SST-c(temp$SST[2:(nrow(temp))],NA)
  temp$sst.diff         <- temp$sst.step.before*temp$sst.step.after
  #temp$wetdry.prob      <- 1-temp$wetdry_ratio/24
  temp$sst.diff2        <- temp$sst.diff#*(temp$wetdry.prob+1)
  temp$SST.remove       <- F
  
  temp$SST.remove[temp$sst.diff2>12] <- T
  
  output <- subset(temp,select=c(dtime,doy,month,year,jday,SST,SST.remove))
  
}