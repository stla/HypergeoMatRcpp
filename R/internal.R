#' @useDynLib HypergeoMatRcpp
#' @importFrom Rcpp evalCpp
#' @importFrom HypergeoMat mvgamma mvbeta
NULL

isPositiveInteger <- function(m){
  is.vector(m) && is.numeric(m) && length(m) == 1L && floor(m) == m
}

isNotNegativeInteger <- function(z){
  Im(z) != 0 || Re(z)>0 || Re(z) != trunc(Re(z))
}

isSymmetricPositive <- function(M){
  isSymmetric(M) && all(eigen(M, symmetric = TRUE, only.values = TRUE)$values >= 0)
}

