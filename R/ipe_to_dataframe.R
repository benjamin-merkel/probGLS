#' transform .ipe file into data frame
#' 
#' transform migrate tech .ipe files into a data frame
#' @param ipefile the ipe file to be transformed
#' @return A data.frame for further use in GeoLight or ProMM
#' @return modified after trnTrans() in GeoLight by Tamara Emmenegger
#' @export


ipe_to_dataframe<-function(ipefile){
  data    <- read.table(ipefile)
  tFirst  <- vector("numeric",length=(length(data[,1])-1))
  tSecond <- vector("numeric",length=(length(data[,1])-1))
  type    <- vector("numeric",length=(length(data[,1])-1))
  ConvInt <- vector("numeric",length=(length(data[,1])-1))
  
  for (i in 1:(length(data[,1])-1)){
    date1      <- as.Date(substr(as.character(data[,1][i]),1,10),format="%Y-%m-%d")
    tFirst[i]  <- as.POSIXct(paste(as.character(date1),
                                   prefix=substr(as.character(data[,2][i]),1,8)),tz="UTC")
    date2      <- as.Date(substr(as.character(data[,1][i+1]),1,10),format="%Y-%m-%d")
    tSecond[i] <- as.POSIXct(paste(as.character(date2),
                                   prefix=substr(as.character(data[,2][i+1]),1,8)),tz="UTC")
    type[i]    <- 0
    
    if(as.character(data[,3][i])=="Sunrise" & as.character(data[,3][i+1])=="Sunset" ) type[i] <- 1
    if(as.character(data[,3][i])=="Sunset"  & as.character(data[,3][i+1])=="Sunrise") type[i] <- 2
    
    ConvInt[i] <- 9
    for(k in 8:1){
      if(data[,5][i+1]==k) ConvInt[i] <- k
      if(data[,5][i]  ==k) ConvInt[i] <- k
    }
  }
  output <- data.frame(tFirst=as.POSIXlt(tFirst,origin="1970-01-01",tz="UTC"),tSecond=as.POSIXlt(tSecond,origin="1970-01-01",tz="UTC"),type=type,ConvInt=ConvInt)
  return(output)
}

