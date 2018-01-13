#' Likelihood functions for negative binomial and mixed poisson-beta
#'
#' The likelihood functions are required for optimisation purposes.
#' There are functions to compute the likelihood of the negative
#' binomial and the mixed poisson-beta distributions with and
#' without a zero inflation.
#'
#'
#' @details
#' Functions Loglik_pois, Loglik_nb, Loglik_mpb compute the negative
#' log-likelihood of poisson, negative binomial and the mixed poisson-beta
#' distributions against the data.
#' Functions Loglik_pois_zero, Loglik_nb_zero, Loglik_mpb_zero compute the
#' zero-inflated log-likelihood values.


#' @param data The dataset for which the likelihood is computed
#' @param par.pois Parameters for the Poisson distribution
#' @param par.nb Parameters for the negative binomial distribution
#' @param par.mpb Parameters for the mixed-poisson-beta distribution
#' @param par.pois.zero,par.nb.zero,par.mpb.zero Parameters for the
#'     respective zero-inflated distributions
#' @param par.pois2,par.nb2,par.mpb2 Parameters for the respective
#'     two population distributions
#'
#' @keywords likelihood negative binomial mixed poisson-beta
#'
#' @name likelihood-nb-mpb
#' @importFrom stats dpois dnbinom rnorm
#' @export
nLoglik_pois <- function(data, par.pois) {
  if (par.pois < 0) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  } else {
    x <- -sum(dpois(x = data, lambda = par.pois, log = TRUE))
    if (x == Inf)
      return (100000 + (rnorm(1, 10000, 20) ^ 2))
    else
      return(x)
  }
}


#' @rdname likelihood-nb-mpb
#' @export
nLoglik_nb <- function(data, par.nb) {
  if (par.nb[1] < 0 || par.nb[2] < 0) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  } else {
    if (sum(log(dnbinom(
      x = data, size = par.nb[1], mu = par.nb[2]
    ))) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else
      return(-sum(log(
        dnbinom(x = data, size = par.nb[1], mu = par.nb[2])
      )))
  }
}


#' @rdname likelihood-nb-mpb
#' @export
nLoglik_mpb <- function(data, par.mpb) {
  if (par.mpb[1] < 0 ||
      par.mpb[2] < 0 || par.mpb[3] < 0) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  } else {
    if (sum(log((
      mpb2::dmpb(
        x = data,
        alpha = par.mpb[1],
        beta = par.mpb[2],
        c = par.mpb[3]
      )
    ))) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else
      return(-sum(log((
        mpb2::dmpb(
          x = data,
          alpha = par.mpb[1],
          beta = par.mpb[2],
          c = par.mpb[3]
        )
      ))))
  }
}

#' @rdname likelihood-nb-mpb
#' @export
nLoglik_pois_zero <- function(data, par.pois.zero) {
  if (par.pois.zero[1] < 0 ||
      par.pois.zero[2] < 0 ||
      par.pois.zero[2] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.pois.zero[2] * (data == 0) + (1 - par.pois.zero[2]) * dpois(x = data, lambda = par.pois.zero[1])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.pois.zero[2] * (data == 0) + (1 - par.pois.zero[2]) * dpois(x = data, lambda = par.pois.zero[1])
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
nLoglik_nb_zero <- function(data, par.nb.zero) {
  if (par.nb.zero[1] < 0 ||
      par.nb.zero[2] < 0 ||
      par.nb.zero[3] < 0 ||
      par.nb.zero[3] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.nb.zero[3] * (data == 0) + (1 - par.nb.zero[3]) * dnbinom(x = data, size = par.nb.zero[1], mu = par.nb.zero[2])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.nb.zero[3] * (data == 0) + (1 - par.nb.zero[3]) * dnbinom(
          x = data,
          size = par.nb.zero[1],
          mu = par.nb.zero[2]
        )
      )))
    }
  }
}


#' @rdname likelihood-nb-mpb
#' @export
nLoglik_mpb_zero <- function(data, par.mpb.zero) {
  if (par.mpb.zero[1] < 0 ||
      par.mpb.zero[2] < 0 ||
      par.mpb.zero[2] < par.mpb.zero[1]  ||
      par.mpb.zero[3] < 0 ||
      par.mpb.zero[4] < 0 ||
      par.mpb.zero[4] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.mpb.zero[4] * (data == 0) + (1 - par.mpb.zero[4]) * mpb2::dmpb(
        x = data,
        alpha = par.mpb.zero[1],
        beta = par.mpb.zero[2],
        c = par.mpb.zero[3]
      )
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.mpb.zero[4] * (data == 0) + (1 - par.mpb.zero[4]) * mpb2::dmpb(
          x = data,
          alpha = par.mpb.zero[1],
          beta = par.mpb.zero[2],
          c = par.mpb.zero[3]
        )
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
nLoglik_pois_two <- function(data, par.pois2) {
  if (par.pois2[1] < 0 ||
      par.pois2[2] < 0 ||
      par.pois2[3] < 0 ||
      par.pois2[3] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.pois2[3] * dpois(x = data, lambda = par.pois2[1]) + (1 - par.pois2[3]) * dpois(x = data, lambda = par.pois2[2])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.pois2[3] * dpois(x = data, lambda = par.pois2[1]) + (1 - par.pois2[3]) * dpois(x = data, lambda = par.pois2[2])
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
nLoglik_nb_two <- function(data, par.nb2) {
  if (par.nb2[1] < 0 ||
      par.nb2[2] < 0 ||
      par.nb2[3] < 0 ||
      par.nb2[4] < 0 ||
      par.nb2[5] < 0 ||
      par.nb2[5] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.nb2[5] * dnbinom(x = data, size = par.nb2[1], mu = par.nb2[2]) + (1 - par.nb2[5]) * dnbinom(x = data, size = par.nb2[3], mu = par.nb2[4])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.nb2[5] * dnbinom(
          x = data,
          size = par.nb2[1],
          mu = par.nb2[2]
        ) + (1 - par.nb2[5]) * dnbinom(
          x = data,
          size = par.nb2[3],
          mu = par.nb2[4]
        )
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
nLoglik_mpb_two <- function(data, par.mpb2) {
  if (par.mpb2[7] < 0 ||
      par.mpb2[2] < 0 ||
      par.mpb2[3] < 0 ||
      par.mpb2[4] < 0 ||
      par.mpb2[5] < 0 ||
      par.mpb2[6] < 0 ||
      par.mpb2[1] < 0 ||
      par.mpb2[1] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.mpb2[1] * mpb2::dmpb(
        x = data,
        alpha = par.mpb2[2],
        beta = par.mpb2[3],
        c = par.mpb2[4]
      ) + (1 - par.mpb2[1]) * mpb2::dmpb(
        x = data,
        alpha = par.mpb2[5],
        beta = par.mpb2[6],
        c = par.mpb2[7]
      )
    )) == -Inf)
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.mpb2[1] * mpb2::dmpb(
          x = data,
          alpha = par.mpb2[2],
          beta = par.mpb2[3],
          c = par.mpb2[4]
        ) + (1 - par.mpb2[1]) * mpb2::dmpb(
          x = data,
          alpha = par.mpb2[5],
          beta = par.mpb2[6],
          c = par.mpb2[7]
        )
      )))
    }
  }
}
