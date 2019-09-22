#' Hypergeometric function of a matrix argument
#'
#' @description Evaluates a truncated hypergeometric function of a matrix
#' argument.
#'
#' @param m truncation weight of the summation, a positive integer
#' @param a the "upper" parameters, a numeric or complex vector,
#' possibly empty (or \code{NULL})
#' @param b the "lower" parameters, a numeric or complex vector,
#' possibly empty (or \code{NULL})
#' @param x either a real or complex square matrix,
#' or a numeric or complex vector, the eigenvalues of the matrix
#' @param alpha the alpha parameter, a positive number
#'
#' @return A real or a complex number.
#' @export
#'
#' @note The hypergeometric function of a matrix argument is usually defined
#' for a symmetric real matrix or a Hermitian complex matrix.
#'
#' @references Plamen Koev and Alan Edelman.
#' \emph{The Efficient Evaluation of the Hypergeometric Function of a Matrix Argument}.
#' Mathematics of Computation, 75, 833-846, 2006.
#'
#' @examples # a scalar x example, the Gauss hypergeometric function
#' hypergeomPFQ(m = 20, a = c(1,2), b = c(3), x = 0.5)
#' gsl::hyperg_2F1(1, 2, 3, 0.5)
#' # 0F0 is the exponential of the trace
#' X <- toeplitz(c(3,2,1))/10
#' hypergeomPFQ(m = 10, a = NULL, b = NULL, x = X)
#' exp(sum(diag(X)))
#' # 1F0 is det(I-X)^(-a)
#' X <- toeplitz(c(3,2,1))/100
#' hypergeomPFQ(m = 15, a = 3, b = NULL, x = X)
#' det(diag(3)-X)^(-3)
#' # Herz's relation for 1F1
#' hypergeomPFQ(m=15, a = 2, b = 3, x = X)
#' exp(sum(diag(X))) * hypergeomPFQ(m=15, a = 3-2, b = 3, x = -X)
#' # Herz's relation for 2F1
#' hypergeomPFQ(15, a = c(1,2), b = 3, x = X)
#' det(diag(3)-X)^(-2) *
#'   hypergeomPFQ(15, a = c(3-1,2), b = 3, -X%*%solve(diag(3)-X))
hypergeomPFQ <- function(m, a, b, x, alpha = 2){
  stopifnot(
    isPositiveInteger(m),
    is.null(a) || (is.vector(a) && is.atomic(a)),
    is.null(b) || (is.vector(b) && is.atomic(b)),
    is.vector(alpha) && is.atomic(alpha),
    length(alpha) == 1L,
    is.numeric(alpha),
    alpha > 0
  )
  if(is.matrix(x)){
    x <- eigen(x, only.values = TRUE)$values
  }else{
    stopifnot(is.atomic(x), is.numeric(x) || is.complex(x))
  }
  if(is.null(a)){
    a <- numeric(0)
  }
  if(is.null(b)){
    b <- numeric(0)
  }
  # if(all(x == x[1L])){
  #   return(HypergeoI(m, alpha, a, b, length(x), x[1L]))
  # }
  #
  Rcpp_hypergeomPFQ(m = as.integer(m), a = a, b = b, x = x, alpha = alpha)
}