% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{load_landmask}
\alias{load_landmask}
\title{Internal function to extract on/off land mask from NOAA OISST V2}
\usage{
load_landmask(FILE_NAME, LONS, LATS)
}
\arguments{
\item{FILE_NAME}{landmask file name with path}

\item{LONS}{Point longitudes}

\item{LATS}{Point latitudes}
}
\description{
This internal function determines if points are on land using the NOAA OISST V2 land maks file.
}
\keyword{internal}
