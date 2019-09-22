#' @useDynLib HypergeoMatRcpp
#' @importFrom Rcpp evalCpp
NULL

isPositiveInteger <- function(m){
  is.vector(m) && is.numeric(m) && length(m) == 1L && floor(m) == m
}

isNotNegativeInteger <- function(z){
  Im(z) != 0 || Re(z)>0 || Re(z) != trunc(Re(z))
}
