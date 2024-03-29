% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trn_to_dataframe.R
\name{trn_to_dataframe}
\alias{trn_to_dataframe}
\title{transform .trn file into data frame}
\usage{
trn_to_dataframe(trnfile)
}
\arguments{
\item{trnfile}{the trn file to be transformed}
}
\value{
A data.frame for further use in GeoLight or probGLS

modified after trnTrans() in GeoLight by Tamara Emmenegger
}
\description{
transform biotrack or BAS .trn files into a data frame with tFirst, tSecond and type
}
\examples{
###########################################
# example black browed albatross light data 
###########################################

# light data ----
trn           <- twilightCalc(BBA_lux$dtime, BBA_lux$lig, ask = FALSE, LightThreshold = 2)

summary(trn)
}
