#' Mixed Poisson-Beta Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the mixed poisson-beta distribution. This is basically a
#' poisson distribution where the parameter itself is again beta distributed
#' with parameters alpha, beta and additionally scaled on (0, c)
#' in contrast to the usual scaling on (0,1).
#'


#' @param x,q  vector of quantiles
#' @param p  vector of probabilities
#' @param n  number of observations
#' @param alpha  first non-negative parameter of the beta distribution (shape1)
#' @param beta  second non-negative parameter of the beta distribution b(shape2)
#' @param c  third parameter of the beta distribution (i.e. the standard
#' beta(alpha,beta) distribution is scaled not on (0, 1), but on (0, c))
#' @param log,log.p  logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail  logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' @keywords mixed poisson-beta distribution

#' @name mpb2
#' @useDynLib mpb2
#' @importFrom Rcpp evalCpp sourceCpp
#' @export
#' @examples
#'  X <- dmpb(x=0:200, alpha=5, beta=3, c=20)
#'  plot(0:200, X, type='l')
#'  Y <- dmpb(0:10, seq(10.0,11.0,by=0.1), seq(30.0,31.0,by=0.1), seq(10.2,11.2,by=0.1))
dmpb <- function(x, alpha, beta, c, log = FALSE) {
  cpp_dmpb(x, alpha, beta, c, log)
}


#' @rdname mpb2
#' @export
#' @examples
#'  Y <- pmpb(q= 0 :200, alpha=5, beta= 3, c=20)
#'  plot(0:200, Y, type="l")
pmpb <- function(q, alpha, beta, c, lower.tail = TRUE, log.p = FALSE) {
  cpp_pmpb(q, alpha, beta, c, lower.tail, log.p)
}


#' @rdname mpb2
#' @export
#' @examples
#'  Z <- qmpb(p= seq(0,1, by= 0.01), alpha=5, beta= 3, c=20)
#'  plot(seq(0,1, by= 0.01),Z, type="l")
qmpb <- function(p, alpha, beta, c, lower.tail = TRUE, log.p = FALSE) {
  cpp_qmpb(p, alpha, beta, c, lower.tail, log.p)
}


#' @rdname mpb2
#' @export
#' @examples
#'  RV <- rmpb(n = 1000, alpha=5, beta= 3, c=20)
#'  plot(0 : 200, X, type="l")
#'  lines(density(RV), col="red")
#'  R2 <- rmpb(11, seq(10.0,11.0,by=0.1), seq(30.0,31.0,by=0.1), seq(10.2,11.2,by=0.1))
rmpb <- function(n, alpha, beta, c) {
  cpp_rmpb(n, alpha, beta, c)
}
