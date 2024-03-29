#' Type one Bessel function of Herz
#'
#' @description Evaluates the type one Bessel function of Herz.
#'
#' @param m truncation weight of the summation, a positive integer
#' @param nu the order parameter, real number with \code{nu > -1}
#' @param x either a real or complex square matrix with real eigenvalues, 
#' or a numeric vector, the eigenvalues of the matrix
#'
#' @return A real number.
#' @export
#'
#' @note This function is usually defined
#' for a symmetric real matrix or a Hermitian complex matrix.
#'
#' @references A. K. Gupta and D. K. Nagar.
#' \emph{Matrix variate distributions}. Chapman and Hall, 1999.
#'
#' @examples # for a scalar x, the relation with the Bessel J-function:
#' t <- 2
#' nu <- 3
#' besselJ(t, nu)
#' BesselA(m=15, t^2/4, nu) * (t/2)^nu
BesselA <- function(m, x, nu){
  stopifnot(
    is.vector(nu) && is.atomic(nu),
    length(nu) == 1L,
    is.numeric(nu),
    Re(nu) > -1
  )
  if(is.matrix(x)){
    p <- nrow(x)
  }else{
    p <- length(x)
  }
  hypergeomPFQ(m, NULL, nu+(p+1)/2, -x) / mvgamma(nu+(p+1)/2, p)
}
