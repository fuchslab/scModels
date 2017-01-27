#' Mixed Poisson-Beta Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the mixed poisson-beta distribution. This is basically a
#' poisson distribution where the parameter itself is again beta distributed
#' with parameters alpha and beta and additional scaled on (0, c).
#' In difference to the usual scaling on (0,1)
#'


#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param alpha first non-negative parameter of the beta distribution (shape1)
#' @param beta second non-negative parameter of the beta distribution b(shape2)
#' @param c thirs parameter of the beta distribution (i.e. the standard
#' beta(alpha,beta) distribution is scaled not on (0, 1), but on (0, c.value))
#' @keywords mixed poisson-beta distribution

#' @name mpb2
#' @useDynLib mpb2
#' @importFrom Rcpp evalCpp
NULL


