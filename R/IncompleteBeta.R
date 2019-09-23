#' Incomplete Beta function of a matrix argument
#'
#' @description Evaluates the incomplete Beta function of a matrix argument.
#'
#' @param m truncation weight of the summation, a positive integer
#' @param a,b real parameters with \code{a > (p-1)/2},
#' \code{b > (p-1)/2}, where \code{p} is the dimension (the order of the matrix)
#' @param x either a real positive symmetric matrix or a complex positive
#' Hermitian matrix "smaller" than the identity matrix (i.e. \code{I-x} is positive),
#' or a numeric vector, the eigenvalues of the matrix
#'
#' @return A real number.
#' @export
#'
#' @references A. K. Gupta and D. K. Nagar.
#' \emph{Matrix variate distributions}. Chapman and Hall, 1999.
#'
#' @examples # for a scalar x, this is the incomplete Beta function:
#' a <- 2; b <- 3
#' x <- 0.75
#' IncBeta(m = 15, a, b, x)
#' gsl::beta_inc(a, b, x)
IncBeta <- function(m, a, b, x){
  if(is.matrix(x)){
    stopifnot(isSymmetricPositive(x))
    p <- nrow(x)
    stopifnot(isSymmetricPositive(diag(p)-x))
  }else if(is.atomic(x)){
    p <- length(x)
  }else{
    stop("Invalid `x` argument")
  }
  stopifnot(
    is.numeric(a),
    a > (p-1)/2,
    is.numeric(b),
    b > (p-1)/2
  )
  if(is.matrix(x)){
    DET <- det(x)
  }else{
    DET <- prod(x)
  }
  DET^a / mvbeta(a, b, p) / a *
    hypergeomPFQ(m, c(a, -b+(p+1)/2), a+(p+1)/2, x)
}
