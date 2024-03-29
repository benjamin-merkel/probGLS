% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sst_deduction.R
\name{sst_deduction}
\alias{sst_deduction}
\title{SST deduction}
\usage{
sst_deduction(datetime, temp, temp.range = c(-2, 30), plot = T)
}
\arguments{
\item{datetime}{date time object as POSIXct in UTC}

\item{temp}{temperature readings}

\item{temp.range}{min and max of temperature range}

\item{plot}{Logical. Plot the results if TRUE}
}
\description{
SST deduction from logger temperature reading
}
\details{
only works with temp data that are recorded while the logger is submerged in seawater
}
\examples{
#################################################
# example black browed albatross temperature data 
#################################################

# sst data ----
sen           <- sst_deduction(datetime = BBA_sst$dtime, temp = BBA_sst$temp, temp.range = c(-2,30))

summary(sen)
}
