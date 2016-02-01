#' transform lotek day log data into data frame
#' 
#' transform lotek daylog file data into a data frame with tFirst, tSecond and type
#' @param date column "TimeS" in daylog file as as.Date()
#' @param sunrise column "Sunrise" in daylog file as as.character()
#' @param sunset column "Sunset" in daylog file as as.character()
#' @param TRLon threshold method calculated longitude in daylog file
#' @param TRLat threshold method calculated latitude in daylog file
#' @return A data.frame for further use in GeoLight or ProMM
#' @return modified after trnTrans() in GeoLight by Tamara Emmenegger
#' @export

lotek_to_dataframe<-function(date,sunrise,sunset,TRLon,TRLat){
  
  data <- data.frame(date,sunrise,sunset,TRLon,TRLat)
  
  data$TRLon[data$TRLon==200] <- NA
  data$TRLat[data$TRLat==100] <- NA
  
  data$sunrise <- as.character(data$sunrise)
  data$sunset  <- as.character(data$sunset)
  
  data$sunrise[data$sunrise==""] <- "12:00"
  data$sunset [data$sunset ==""] <- "12:00"
  
  data$sr.H <- as.numeric(matrix(unlist(strsplit(as.character(data$sunrise),':')),ncol=2,byrow=T)[,1])
  data$sr.M <- as.numeric(matrix(unlist(strsplit(as.character(data$sunrise),':')),ncol=2,byrow=T)[,2])
  
  data$ss.H <- as.numeric(matrix(unlist(strsplit(as.character(data$sunset),':')),ncol=2,byrow=T)[,1])
  data$ss.M <- as.numeric(matrix(unlist(strsplit(as.character(data$sunset),':')),ncol=2,byrow=T)[,2])
  
  data$sr <- data$date
  data$ss <- data$date
  
  data$sr[data$sr.H>23] <- data$sr[data$sr.H>23]+1
  data$ss[data$ss.H>23] <- data$ss[data$ss.H>23]+1
  
  data$sr.H[data$sr.H>23] <- data$sr.H[data$sr.H>23]-24
  data$ss.H[data$ss.H>23] <- data$ss.H[data$ss.H>23]-24
  
  data$sunrise2 <- as.POSIXct(strptime(paste(data$sr,data$sr.H,data$sr.M,sep=" "),'%Y-%m-%d %H %M'),tz="UTC")
  data$sunset2  <- as.POSIXct(strptime(paste(data$ss,data$ss.H,data$ss.M,sep=" "),'%Y-%m-%d %H %M'),tz="UTC")
  
  times         <- data.frame(tFirst = data$sunrise2)
  times$tSecond <- data$sunset2
  times$tFirst  <- twilight(tm = times$tFirst , lon = data$TRLon, lat = data$TRLat, rise = TRUE , zenith = 93.44,iter=5) 
  times$tSecond <- twilight(tm = times$tSecond, lon = data$TRLon, lat = data$TRLat, rise = FALSE, zenith = 93.44,iter=5)
  
  tl                               <- nrow(times)
  data2                            <- rbind(times,times)
  data2$V2                         <- c(rep('Sunrise',length=tl),rep('Sunset',length=tl))
  data2$tFirst[data2$V2=="Sunset"] <- data2$tSecond[data2$V2=="Sunset"]
  data2                            <- subset(data2,select=c(tFirst,V2))
  colnames(data2)                  <- c("V1","V2")
  data2                            <- data2[order(data2$V1),]
  data2$dtime                      <- as.Date(data2$V1)
  
  tFirst  <- vector("numeric",length=(nrow(data2)-1))
  tSecond <- vector("numeric",length=(nrow(data2)-1))
  type    <- vector("numeric",length=(nrow(data2)-1))
  
  for (idata in 1:(nrow(data2)-1)){
    tFirst[idata]  <- data2$V1[idata]
    tSecond[idata] <- data2$V1[idata+1]
    type[idata]    <- 0
    if(as.character(data2$V2[idata])=="Sunrise" & as.character(data2$V2[idata+1])=="Sunset" ) type[idata] <- 1
    if(as.character(data2$V2[idata])=="Sunset"  & as.character(data2$V2[idata+1])=="Sunrise") type[idata] <- 2
  }
  
  output <- data.frame(tFirst  = as.POSIXlt(tFirst ,origin="1970-01-01",tz="UTC"),
                       tSecond = as.POSIXlt(tSecond,origin="1970-01-01",tz="UTC"),
                       type)
  output <- output[output$type>0,]
  
  return(output)
}
