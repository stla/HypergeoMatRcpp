#' Hypergeometric function of a matrix argument
#'
#' @description Evaluates a truncated hypergeometric function of a matrix
#' argument.
#'
#' @param m truncation weight of the summation, a positive integer
#' @param a the "upper" parameters, a numeric vector,
#' possibly empty (or \code{NULL})
#' @param b the "lower" parameters, a numeric vector,
#' possibly empty (or \code{NULL})
#' @param x either a real or complex square matrix with real eigenvalues,
#' or a numeric vector, the eigenvalues of the matrix
#' @param alpha the alpha parameter, a positive number
#'
#' @return A real number.
#' @export
#' 
#' @details This is an implementation of Koev & Edelman's algorithm
#' (see the reference). This algorithm is split into two parts: the case of
#' a scalar matrix (multiple of an identity matrix) and the general case.
#' The case of a scalar matrix is much faster (try e.g. \code{x = c(1,1,1)} vs
#' \code{x = c(1,1,0.999)}).
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
    is.null(a) || (is.vector(a) && is.atomic(a) && is.numeric(a)),
    is.null(b) || (is.vector(b) && is.atomic(b) && is.numeric(b)),
    is.vector(alpha) && is.atomic(alpha),
    length(alpha) == 1L,
    is.numeric(alpha),
    alpha > 0
  )
  if(is.matrix(x)){
    x <- eigen(x, only.values = TRUE)$values
    if(any(is.complex(x))){
      stop("The eigenvalues of `x` are not all real")
    }
  }else{
    stopifnot(is.atomic(x), is.numeric(x))
  }
  if(is.null(a)){
    a <- numeric(0)
  }
  if(is.null(b)){
    b <- numeric(0)
  }
  Rcpp_hypergeomPFQ(m = as.integer(m), a = a, b = b, x = x, alpha = alpha)
}