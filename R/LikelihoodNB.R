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
#' @param par.nb Parameters for the negative binomial distribution
#' @param par.mpb Parameters for the mixed-poisson-beta distribution
#'
#' @keywords likelihood negative binomial mixed poisson-beta
#'
#' @name likelihood-nb-mpb
#' @importFrom stats dnbinom rnorm
#' @export
Loglik_pois <- function(data, par.pois) {
  if (par.pois < 0) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  } else {
    if (sum(log(dpois(x = data, lambda = par.pois))) == -Inf)
      return (100000 + (rnorm(1, 10000, 20) ^ 2))
    else
      return (-sum(log(dpois(
        x = data, lambda = data
      ))))
  }
}


#' @rdname likelihood-nb-mpb
#' @export
Loglik_nb <- function(data, par.nb) {
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
Loglik_mpb <- function(data, par.mpb) {
  if (par.mpb[1] < 0 ||
      par.mpb[2] < 0 || par.mpb[3] < 0 || par.mpb[2] < par.mpb[1]) {
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
Loglik_pois_zero <- function(data, par1_z) {
  if (par1_z[1] < 0 ||
      par1_z[2] < 0 ||
      par1_z[2] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par1_z[2] * (data == 0) + (1 - par1_z[2]) * dpois(x = data, lambda = par1_z[1])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par1_z[2] * (data == 0) + (1 - par1_z[2]) * dpois(x = data, lambda = par1_z[1])
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
Loglik_nb_zero <- function(data, par.nb) {
  if (par.nb[1] < 0 ||
      par.nb[2] < 0 ||
      par.nb[3] < 0 ||
      par.nb[3] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.nb[3] * (data == 0) + (1 - par.nb[3]) * dnbinom(x = data, size = par.nb[1], mu = par.nb[2])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.nb[3] * (data == 0) + (1 - par.nb[3]) * dnbinom(
          x = data,
          size = par.nb[1],
          mu = par.nb[2]
        )
      )))
    }
  }
}


#' @rdname likelihood-nb-mpb
#' @export
Loglik_mpb_zero <- function(data, par.mpb) {
  if (par.mpb[1] < 0 ||
      par.mpb[2] < 0 ||
      par.mpb[2] < par.mpb[1]  ||
      par.mpb[3] < 0 ||
      par.mpb[4] < 0 ||
      par.mpb[4] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par.mpb[4] * (data == 0) + (1 - par.mpb[4]) * mpb2::dmpb(
        x = data,
        alpha = par.mpb[1],
        beta = par.mpb[2],
        c = par.mpb[3]
      )
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par.mpb[4] * (data == 0) + (1 - par.mpb[4]) * mpb2::dmpb(
          x = data,
          alpha = par.mpb[1],
          beta = par.mpb[2],
          c = par.mpb[3]
        )
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
Loglik_pois_two <- function(data, par1_2) {
  if (par1_2[1] < 0 ||
      par1_2[2] < 0 ||
      par1_2[3] < 0 ||
      par1_2[3] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par1_2[3] * dpois(x = data, lambda = par1_2[1]) + (1 - par1_2[3]) * dpois(x = data, lambda = par1_2[2])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par1_2[3] * dpois(x = data, lambda = par1_2[1]) + (1 - par1_2[3]) * dpois(x = data, lambda = par1_2[2])
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
Loglik_nb_two <- function(data, par2_2) {
  if (par2_2[1] < 0 ||
      par2_2[2] < 0 ||
      par2_2[3] < 0 ||
      par2_2[4] < 0 ||
      par2_2[5] < 0 ||
      par2_2[5] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par2_2[5] * dnbinom(x = data, size = par2_2[1], mu = par2_2[2]) + (1 - par2_2[5]) * dnbinom(x = data, size = par2_2[3], mu = par2_2[4])
    )) == -Inf)
      return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par2_2[5] * dnbinom(
          x = data,
          size = par2_2[1],
          mu = par2_2[2]
        ) + (1 - par2_2[5]) * dnbinom(
          x = data,
          size = par2_2[3],
          mu = par2_2[4]
        )
      )))
    }
  }
}

#' @rdname likelihood-nb-mpb
#' @export
Loglik_mpb_two <- function(data, par3_2) {
  if (par3_2[7] < 0 ||
      par3_2[2] < 0 ||
      par3_2[3] < 0 ||
      par3_2[4] < 0 ||
      par3_2[5] < 0 ||
      par3_2[6] < 0 ||
      par3_2[1] < 0 ||
      par3_2[1] > 1) {
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
  }
  else {
    if (sum(log(
      par3_2[1] * mpb2::dmpb(
        x = data,
        alpha = par3_2[2],
        beta = par3_2[3],
        c = par3_2[4]
      ) + (1 - par3_2[1]) * mpb2::dmpb(
        x = data,
        alpha = par3_2[5],
        beta = par3_2[6],
        c = par3_2[7]
      )
    )) == -Inf)
    return(100000 + (rnorm(1, 10000, 20) ^ 2))
    else{
      return(-sum(log(
        par3_2[1] * mpb2::dmpb(
          x = data,
          alpha = par3_2[2],
          beta = par3_2[3],
          c = par3_2[4]
        ) + (1 - par3_2[1]) * mpb2::dmpb(
          x = data,
          alpha = par3_2[5],
          beta = par3_2[6],
          c = par3_2[7]
        )
      )))
    }
  }
}
