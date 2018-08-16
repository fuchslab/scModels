#' Poisson-beta Distribution
#'
#' Density, distribution function, quantile function and random generation for the
#' Poisson-beta distribution: a Poisson distribution whose parameter itself follows
#' a beta distribution. Alpha and beta are the parameters of this specific beta
#' distribution which is scaled on (0, c) in contrast to the usual scaling on (0,1).
#'


#' @param x,q  vector of (non-negative integer) quantiles
#' @param p  vector of probabilities
#' @param n  number of observations
#' @param alpha,beta  non-negative parameters of the beta distribution (shape1 and shape2)
#' @param c  numeric scaling parameter of the beta distribution. The standard beta is
#'     scaled on (0,1) (default) and can be transformed to (0,c).
#' @param log,log.p  logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail  logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' @keywords poisson-beta distribution
#' @name poisson-beta
#' @useDynLib scModels
#' @importFrom Rcpp evalCpp sourceCpp
#' @export
#' @examples
#'  X <- dpb(x=0:200, alpha=5, beta=3, c=20)
#'  plot(0:200, X, type='l')
#'  Y <- dpb(0:10, seq(10.0,11.0,by=0.1), seq(30.0,31.0,by=0.1), seq(10.2,11.2,by=0.1))
dpb <- function(x, alpha, beta, c = 1, log = FALSE) {
  cpp_dpb(x, alpha, beta, c, log)
}


#' @rdname poisson-beta
#' @export
#' @examples
#'  Y <- ppb(q= 0 :200, alpha=5, beta= 3, c=20)
#'  plot(0:200, Y, type="l")
ppb <- function(q, alpha, beta, c = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_ppb(q, alpha, beta, c, lower.tail, log.p)
}


#' @rdname poisson-beta
#' @export
#' @examples
#'  Z <- qpb(p= seq(0,1, by= 0.01), alpha=5, beta= 3, c=20)
#'  plot(seq(0,1, by= 0.01),Z, type="l")
qpb <- function(p, alpha, beta, c = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qpb(p, alpha, beta, c, lower.tail, log.p)
}


#' @rdname poisson-beta
#' @export
#' @examples
#'  RV <- rpb(n = 1000, alpha=5, beta= 3, c=20)
#'  plot(0 : 200, X, type="l")
#'  lines(density(RV), col="red")
#'  R2 <- rpb(11, seq(10.0,11.0,by=0.1), seq(30.0,31.0,by=0.1), seq(10.2,11.2,by=0.1))
rpb <- function(n, alpha, beta, c = 1) {
  cpp_rpb(n, alpha, beta, c)
}
