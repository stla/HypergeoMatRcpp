#' @useDynLib HypergeoMatRcpp
#' @importFrom Rcpp evalCpp
#' @importFrom HypergeoMat mvgamma mvbeta
NULL

isPositiveInteger <- function(m){
  is.vector(m) && is.numeric(m) && length(m) == 1L && floor(m) == m
}

isSymmetricPositive <- function(M){
  isSymmetric(M) && all(eigen(M, symmetric = TRUE, only.values = TRUE)$values >= 0)
}

