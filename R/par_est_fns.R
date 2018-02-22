#' Functions to compute parameter estimates
#'
#' @param x Gene expression data as an array
#' @param iter number of bootstrap replicates to estimate initial
#'     parameters for the mpb
#' @param type keyword for the distribution the data is to be fitted
#'     against. Possible values are ("pois", "zip", "nb", "zinb",
#'     "mpb", "zimpb", "nb2", "mpb2")
#' @keywords parameter estimation
#' @name par-est-fns
#' @importFrom stats kmeans optim
#' @export
estimate_mpb_optim_init <- function(x, iter = 200) {
  sampled_params <- c()
  n <- length(x)
  for (i in 1:iter) {
    d1 <- sample(x, n, replace = TRUE)
    e1 <- mean(d1)
    e2 <- mean(d1 * ( d1 - 1 ))
    e3 <- mean(d1 * (d1 - 1) * (d1 - 2))
    r1 <- e1
    r2 <- e2 / e1
    r3 <- e3 / e2
    x1 <- r1 * r2 - 2 * r1 * r3 + r2 * r3
    x2 <- r1 - 2 * r2 + r3
    alpha <- 2 * r1 * (r3 - r2) / x1
    if (alpha > 0)
      sampled_params <- c(sampled_params, alpha)
  }
  cm <- c(mean(sampled_params), 0, max(data))
  cm[2] <- (function(a, c, m) a * c / m - a)(cm[1], cm[3], mean(x))
  return(cm)
}


#' @rdname par-est-fns
#' @export
get_0inf_parameter <- function(x) length(c(which(x == 0))) / length(x)


#' @rdname par-est-fns
#' @export
get_fitted_params <- function(x, type) {
  if (type == "pois") {
    p <- c(0)
    t <- system.time(o <- optim(par = p, fn = nLoglik_pois, data = x))
  } else if (type == "zip") {
    p <- c(get_0inf_parameter(x), 1)
    t <- system.time(o <- optim(par = p, fn = nLoglik_pois_zero, data = x))
  } else if (type == "nb") {
    p <- c(1, 1)
    t <- system.time(o <- optim(par = p, fn = nLoglik_nb, data = x))
  } else if (type == "zinb") {
    p <- c(get_0inf_parameter(x), 1, 1)
    t <- system.time(o <- optim(par = p, fn = nLoglik_nb_zero, data = x))
  } else if (type == "mpb") {
    p <- estimate_mpb_optim_init(x)
    t <- system.time(o <- optim(par = p, fn = nLoglik_mpb, data = x, control = list(reltol = 0.001, maxit = 100)))
  } else if (type == "zimpb") {
    p <- c(get_0inf_parameter(x), estimate_mpb_optim_init(x))
    t <- system.time(o <- optim(par = p, fn = nLoglik_mpb_zero, data = x, control = list(reltol = 0.001, maxit = 200)))
  } else if (type == "pois2") {
    k <- kmeans(x, centers = 2)
    p <- c(length(which(k$cluster == 1))/length(x), 1, 1)
    t <- system.time(o <- optim(par = p, fn = nLoglik_pois_two, data = x))
  } else if (type == "nb2") {
    k <- kmeans(x, centers = 2)
    p <- c(length(which(k$cluster == 1))/length(x), 1, 1, 1, 1)
    t <- system.time(o <- optim(par = p, fn = nLoglik_nb_two, data = x))
  } else if (type == "mpb2") {
    k <- kmeans(x = x, centers = 2)
    c1 <- x[which(k$cluster == 1)]
    c2 <- x[which(k$cluster == 2)]
    t1 <- estimate_mpb_optim_init(c1)
    t2 <- estimate_mpb_optim_init(c2)
    p <- length(c1)/length(x)
    par <- c(p,t1, t2)
    t <- system.time(o <- optim(par = par, fn = nLoglik_mpb_two, data = x, control = list(reltol = 0.001, maxit = 500)))
  } else {
    warning("Invalid distribution type.")
    return(NULL)
  }
  fit_param <- o
  fit_param$time <- t
  fit_param$AIC <- 2 * length(o$par) + 2 * o$value
  fit_param$BIC <- log(length(x)) * length(o$par) + 2 * o$value
  return(fit_param)
}
