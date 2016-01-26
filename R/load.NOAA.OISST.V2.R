#' load NOAA OISST V2
#' 
#' This function takes 1-year-long NetCDF files of daily SST from the
#' ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/ directory
#' where filenames for daily mean SST files are named with the scheme
#' sst.day.mean.YEAR.v2.nc
#' @param fname full path to NetCDF data file
#' @param lsmask full path to land-sea mask NetCDF file
#' @param lonW western-most longitude of search area, must be smaller than lonE
#' @param lonE eastern-most longitude of search area, must be larger than lon1
#' @param latS southern-most latitude of search area, must be smaller than latN
#' @param latN northern-most latitude of search area, must be larger than latS
#' @param date1 first date in file to extract, must be Date class
#' @param date2 last date in file to extract, must be Date class
#' @param use.landmask use land mask
#' @param extract.value which data to extract: "sst" - SST, "err" - SST error, "icec" - sea ice concentration
#' @return A 3-dimensional array with latitudes in rows, longitudes in columns, and dates along the 3rd dimension. The value [1,1,1] is the northernmost, westernmost lat/long location on the 1st date. The value [1,1,2] is the 2nd date at the same lat/long location (if more than 1 date is requested).
#' @return To extract lat/lon/date values from the output array, use the dimnames() function:
#' @return lats = as.numeric(dimnames(sst2)$Lat)
#' @return longs = as.numeric(dimnames(sst2)$Long)
#' @return dates = as.Date(dimnames(sst2)$Date)
#' @return  
#' @return NetCDF files should be downloaded from the links on:
#' @return http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
#' @return In addition to the temperature data files, also download a copy of the landmask file lsmask.oisst.v2.nc from the same page. Inside the NetCDF files, data are available on a 0.25 degree latitude x 0.25 degree longitude global grid (720x1440 cells) From -89.875N to 89.875N, 0.125E to 359.875E. Locations are at the CENTER of a grid cell.
#' @return modified after Luke Miller Nov 25, 2014; https://github.com/millerlp/Misc_R_scripts/blob/master/NOAA_OISST_ncdf4.R
#' @export



# taken from the script FUNCTION_load.sst

load.NOAA.OISST.V2 = function(fname,lsmask,lonW,lonE,latS,latN, 
                             date1, date2,use.landmask=F,
                             extract.value='sst'){
  
  # Generate set of grid cell latitudes (center of cell) from south to north
  lats = seq(-89.875,89.875,0.25)
  # Generate set of grid cell longitudes (center of cell)
  lons = seq(0.125,359.875,0.25)
  # Create connection to NetCDF data file
  nc = nc_open(fname)
  lonWindx = which.min(abs(lonW - lons)) #get index of nearest longitude value
  if (missing(lonE)){
    # If lonE isn't specified, reused lon1
    lonE = lonW
    lonEindx = lon1indx
    cat("Only 1 longitude specified\n")
  } else {
    # Get index of nearest longitude value to lonE
    lonEindx = which.min(abs(lonE - lons))
  }
  latSindx = which.min(abs(latS - lats)) #get index of nearest latitude value
  if (missing(latN)){
    # If latN is not specified, reuse latS
    latN = latS
    latNindx = latSindx
    cat("Only 1 latitude specified\n")
  } else {
    # Get index of nearest latitude value to latN
    latNindx = which.min(abs(latN - lats))
  }
  # The lon/lat indx values should now correspond to indices in the NetCDF
  # file for the desired grid cell.
  nlon = (lonEindx - lonWindx) + 1 # get number of longitudes to extract
  nlat = (latNindx - latSindx) + 1 # get number of latitudes to extract
  # Extract available dates from netCDF file
  ncdates = nc$dim$time$vals
  ncdates = as.Date(ncdates,origin = '1800-1-1') #available time points in nc
  if (class(date1) == 'Date'){
    # Get index of nearest time point
    date1indx = which.min(abs(date1 - ncdates))
  } else if (class(date1) == 'character'){
    # Convert to a Date object first
    date1 = as.Date(date1)
    date1indx = which.min(abs(date1 - ncdates))
  }
  if (missing(date2)) {
    # If date2 isn't specified, reuse date1
    date2indx = which.min(abs(date1 - ncdates))
    cat('Only 1 date specified\n')
  } else {
    if (class(date2) == 'Date'){
      # If date2 exists, get index of nearest time point to date2
      date2indx = which.min(abs(date2 - ncdates))
    } else if (class(date2) == 'character'){
      date2 = as.Date(date2)
      date2indx = which.min(abs(date2 - ncdates))
    }
  }
  ndates = (date2indx - date1indx) + 1 #get number of time steps to extract
  # Define the output array
  sstout = matrix(data = NA, ncol = nlon, nrow = nlat)
  # Extract the data from the NetCDF file
  sstout[,] = ncvar_get(nc, varid = extract.value,
                         start = c(lonWindx,latSindx),
                         count = c(nlon,nlat))
  
  # close SST ncdf
  nc_close(nc)
  
  
  # The output array sstout will be arranged with longitudes in rows,
  # increasing in an easterly direction as you move down a row (larger
  # longitude values), and latitudes in columns, increasing in latitude (more
  # northerly) as you move across columns. The 3rd dimension represents
  # different dates. This arrangement stems from how the data are set up in
  # the NetCDF file, so my apologies if it's counterintuitive.
  # If there are missing data in the NetCDF, they should appear as 32767.
  # Replace that value with NA if it occurs anywhere.
  sstout = ifelse(sstout == 32767, NA, sstout)
  # The NOAA OISST files contain sea surface temperatures for the entire
  # globe, including on the continents. This clearly isn't right, so they also
  # supply a land-sea mask file in netCDF format. We use the values (0 or 1)
  # stored in the mask file to turn all of the continent areas into NA's.
  # Open the land-sea mask
  if(use.landmask==T) {
    
    nc2 = nc_open(lsmask)
    # Create array to hold land-sea mask
    mask = matrix(data = NA, dim = c(nlon,nlat))
    # Get land-sea mask values (0 or 1)
    mask[,] = ncvar_get(nc2, varid = "lsmask",
                         start = c(lonWindx,latSindx), count = c(nlon,nlat))
    
    #close land mask
    nc_close(nc2)
    
    
    # Replace 0's with NA's
    mask = ifelse(mask == 0,NA,1)
    # Get dimensions of sstout array
    dims = dim(sstout)
    
    sstout[,] = sstout[,] * mask[,] # All masked values become NA
    
  }
  
  # Add dimension names
  attr(sstout,'dimnames') = list(Long = seq(lons[lonWindx],lons[lonEindx],
                                            by = 0.25),
                                 Lat = seq(lats[latSindx],lats[latNindx],
                                           by = 0.25))
  
  # sstout now has dimension names that show the longitude and latitude of
  # each point in the array, as well as the date (3rd dimension of the array).
  ############################################################################
  # Rearrange the output matrix or array so that latitudes run from north to
  # south down the rows, and longitudes run from west to east across columns.
  dims = dim(sstout) # get size of array
  # Make new output array to hold rearranged data. The dimension names will
  # match the newly rearranged latitude and longitude values
  sst2 = matrix(data = NA, dim = c(dims[2],dims[1]),
               dimnames = list(Lat = rev(seq(lats[latSindx],lats[latNindx],
                                             by = 0.25)),
                               Long = seq(lons[lonWindx],lons[lonEindx],by = 0.25)))
  # Step through each page of array and rearrange lat/lon values
   # Extract one day's worth of lat/lon pairs
    temp = as.matrix(sstout[,])
    temp = t(temp) # transpose lon/lat to lat/lon
    temp = temp[nrow(temp):1,] # reverse row order to reverse latitudes
    sst2[,] = temp # write data to sst2 array
  
  ##########################
  sst2 # return sst2 array
  ##########################
  
  
} # end of function
