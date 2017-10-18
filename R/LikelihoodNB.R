#' Likelihood functions for negative binomial and mixed poisson-beta
#'
#' The likelihood functions are required for optimisation purposes.
#' There are functions to compute the likelihood of the negative
#' binomial and the mixed poisson-beta distributions with and
#' without a zero inflation.
#'


#' @param data The dataset for which the likelihood is computed
#' @param par.nb Parameters for the negative binomial distribution
#' @param par.mpb Parameters for the mixed-poisson-beta distribution
#'
#' @keywords likelihood negative binomial mixed poisson-beta
#'
#' @name likelihood-nb-mpb
#' @importFrom stats dnbinom rnorm
#' @export
Loglik_nb <- function(data, par.nb) {
    if (par.nb[1] < 0 || par.nb[2] < 0 ) {return(100000 + (rnorm(1,10000, 20)^2))}
    else { if(sum(log(dnbinom(x = data, size = par.nb[1], mu = par.nb[2]))) == -Inf) return(100000 + (rnorm(1,10000, 20)^2))
    else{    return(- sum(log(dnbinom(x = data, size = par.nb[1], mu = par.nb[2]))))
        }
    }
}


#' @rdname likelihood-nb-mpb
#' @export
Loglik_mpb <- function(data, par.mpb) {
  if (par.mpb[1] < 0 || par.mpb[2] < 0 || par.mpb[3] < 0 || par.mpb[2] < par.mpb[1]  ) {return(100000 + (rnorm(1,10000, 20)^2))}
  else {if( sum(log((mpb2::dmpb(x = data, alpha = par.mpb[1], beta = par.mpb[2], c = par.mpb[3])))) == -Inf) return(100000 + (rnorm(1,10000, 20)^2))
    else{return(- sum(log((mpb2::dmpb(x = data, alpha = par.mpb[1], beta = par.mpb[2], c = par.mpb[3])))))
    }
  }
}


#' @rdname likelihood-nb-mpb
#' @export
Loglik_nb_zero <- function(data, par.nb) {
  if (par.nb[1] < 0 || par.nb[2] < 0 || par.nb[3] < 0 || par.nb[3] > 1) {return(100000 + (rnorm(1,10000, 20)^2))}
  else { if(sum(log(par.nb[3] * (data == 0) + (1 - par.nb[3]) * dnbinom(x = data, size = par.nb[1], mu = par.nb[2]))) == -Inf) return(100000 + (rnorm(1,10000, 20)^2))
    else{    return(- sum(log(par.nb[3] * (data == 0) + (1 - par.nb[3]) * dnbinom(x = data, size = par.nb[1], mu = par.nb[2]))))
    }
  }
}


#' @rdname likelihood-nb-mpb
#' @export
Loglik_mpb_zero <- function(data, par.mpb) {
  if (par.mpb[1] < 0 || par.mpb[2] < 0 || par.mpb[2] < par.mpb[1]  || par.mpb[3] < 0 ||par.mpb[4] < 0 || par.mpb[4] > 1) {return(100000 + (rnorm(1,10000, 20)^2))}
  else { if(sum(log(par.mpb[4] * (data == 0) + (1 - par.mpb[4]) * mpb2::dmpb(x = data, alpha = par.mpb[1], beta = par.mpb[2], c = par.mpb[3]))) == -Inf) return(100000 + (rnorm(1,10000, 20)^2))
    else{
      return(- sum(log(par.mpb[4] * (data == 0) + (1 - par.mpb[4]) * mpb2::dmpb(x = data, alpha = par.mpb[1], beta = par.mpb[2], c = par.mpb[3]))))
    }
  }
}
