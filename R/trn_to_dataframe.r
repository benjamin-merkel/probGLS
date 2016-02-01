#' transform .trn file into data frame
#' 
#' transform biotrack or BAS .trn files into a data frame with tFirst, tSecond and type
#' @param trnfile the trn file to be transformed
#' @return A data.frame for further use in GeoLight or ProMM
#' @return modified after trnTrans() in GeoLight by Tamara Emmenegger
#' @export

trn_to_dataframe<-function(trnfile){
  data    <- read.table(trnfile,sep=",")
  tFirst  <- vector("numeric",length=(length(data[,1])-1))
  tSecond <- vector("numeric",length=(length(data[,1])-1))
  type    <- vector("numeric",length=(length(data[,1])-1))
  ConvInt <- vector("numeric",length=(length(data[,1])-1))
  
  for (i in 1:(length(data[,1])-1)){
    date1      <- as.Date(substr(as.character(data$V1[i]),1,8),format="%d/%m/%y")
    tFirst[i]  <- as.POSIXct(paste(as.character(date1),
                                  prefix=substr(as.character(data$V1[i]),10,17)),tz="UTC")
    date2      <- as.Date(substr(as.character(data$V1[i+1]),1,8),format="%d/%m/%y")
    tSecond[i] <- as.POSIXct(paste(as.character(date2),
                                   prefix=substr(as.character(data$V1[i+1]),10,17)),tz="UTC")
    type[i]    <- 0
    
    if(as.character(data$V2[i])=="Sunrise" & as.character(data$V2[i+1])=="Sunset" ) type[i] <- 1
    if(as.character(data$V2[i])=="Sunset"  & as.character(data$V2[i+1])=="Sunrise") type[i] <- 2
    
    ConvInt[i] <- 9
    for(k in 8:1){
      if(data$V3[i+1]==k) ConvInt[i] <- k
      if(data$V3[i]  ==k) ConvInt[i] <- k
    }
  }
  output <- data.frame(tFirst=as.POSIXlt(tFirst,origin="1970-01-01",tz="UTC"),tSecond=as.POSIXlt(tSecond,origin="1970-01-01",tz="UTC"),type=type,ConvInt=ConvInt)
  return(output)
}

