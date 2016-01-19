#' transform list into data frame
#' 
#' @param x list to be tranformed
#' @return A data.frame for further use in ProMM

f = function(x) function(i) sapply(x, `[[`, i)
