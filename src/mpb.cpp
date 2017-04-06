// [[Rcpp::depends(RcppGSL)]]

#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include "mpb.h"
#include "shared.h"
using namespace Rcpp;

// kummer series using GSL
double kummer_(double x, double a, double b, bool log_v) {
  if(!validKummerParameters(a, b)) {
    return R_NaN;
  }
  double res = gsl_sf_hyperg_1F1(a, b, x);
  if(log_v) {
    return log(res);
  } else {
    return res;
  }
}

// density function
double dmpb_(double x, double alpha, double beta, double c) {
  if( isInadmissible(x) || isInadmissible(alpha) || isInadmissible(beta) || isInadmissible(c) )
    return x+alpha+beta+c;

  if( !isInteger(x) || x < 0  || traits::is_infinite<REALSXP>(x) )
    return 0;

  double cre = kummer_(-c, alpha+x, beta+alpha+x, 1);
  if(x <= 0) {
    return exp(cre);
  } else {
    int sign = (x-1 < 0) ? -1 : 1;
    int x2 = (x-1 > 0) ? (int)std::floor(x-1) : (int)std::floor(1-x);
    double num = 0, denom = 0;
    for(int i=0; i <= x2; i++) {
      num += log((alpha + sign*i));
      denom += log(alpha + beta + sign * i);
    }
    num += x * log(c);
    denom += lgamma(x+1);
    return exp(num-denom+cre);
  }
}

// distribution function
double pmpb_(double x, double alpha, double beta, double c) {
  if( isInadmissible(x) || isInadmissible(alpha) || isInadmissible(beta) || isInadmissible(c) )
    return x+alpha+beta+c;

  if( !isInteger(x) )
    return 0;
  if(traits::is_infinite<REALSXP>(x))
    return 1;
  double res = 0;
  for(int i = 0; i <= x; i++) {
    res += dmpb_(i, alpha, beta, c);
  }
  return res;
}

// distribution function array
double* pmpb_(double alpha, double beta, double c) {
  double *res = (double *)std::malloc(Q_LIMIT * sizeof(double));
  res[0] = dmpb_(0, alpha, beta, c);
  for(int i = 1; i < Q_LIMIT; i++) {
    res[i] = res[i-1] + dmpb_(i, alpha, beta, c);
  }
  return res;
}

// quantiles for single parameters
double qmpb_(double p, double *p_distr) {
  if(isInadmissible(p))
    return NA_REAL;
  if(!validProbability(p)){
    warning("NaNs produced");
    return R_NaN;
  }

  int i;
  if(p > p_distr[Q_LIMIT-1]) {
    return Q_LIMIT;
  }
  if(p < p_distr[0]) {
    return 0;
  }
  for(i = 1; i < Q_LIMIT; i++) {
    if(p > p_distr[i-1] && p < p_distr[i]) {
      return i;
    }
  }
  return i;
}

// quantiles for vectorised parameters
double qmpb_(double p, double alpha, double beta, double c) {
  if(isInadmissible(p))
    return NA_REAL;
  if(!validProbability(p)){
    warning("NaNs produced");
    return R_NaN;
  }

  double *p_distr = pmpb_(alpha, beta, c);
  int i;
  if(p > p_distr[Q_LIMIT-1]) {
    return Q_LIMIT;
  }
  if(p < p_distr[0]) {
    return 0;
  }
  for(i = 1; i < Q_LIMIT; i++) {
    if(p > p_distr[i-1] && p < p_distr[i]) {
      return i;
    }
  }
  return i;
}

//' Kummer's (confluent hypergeometric) function
//'
//' Kummer's (confluent hypergeometric) function of the first kind
//' for numeric (non-complex) values and input parameters
//' @param x numeric value or vector
//' @param a,b numeric parameters of the Kummer function
//' @param log_v logical; if TRUE, the log of the value is returned
//' @name kummer
//' @rdname kummer
//' @export
// [[Rcpp::export]]
NumericVector chf_1F1_gsl(NumericVector x, NumericVector a, NumericVector b, const bool& log_v = false) {
    if(min(NumericVector::create(x.length(), a.length(), b.length())) < 1) {
      return NumericVector(0);
    }
    int n = max(NumericVector::create(x.length(), a.length(), b.length()));
    NumericVector res(n);
    for(int i = 0; i < n; i++) {
        res[i] = kummer_(GETV(x, i), GETV(a, i), GETV(b, i), log_v);
    }
    return res;
}


// [[Rcpp::export]]
NumericVector cpp_dmpb(NumericVector x, NumericVector alpha, NumericVector beta, NumericVector c, const bool& log_p = false) {
  int n = x.size(), type = INPUT_SINGLE;
  if(1 == alpha.size() && 1 == beta.size() && 1 == c.size()) {
    type = INPUT_SINGLE;
  } else if(n == alpha.size() && n == beta.size() && n == c.size()) {
    type = INPUT_VECTORISED;
  } else {
    warning("Dimensions do not match");
  }
  NumericVector res(n);
  for(int i = 0; i < n; i++) {
    switch(type) {
    case INPUT_SINGLE:
      res[i] = dmpb_(x[i], alpha[0], beta[0], c[0]);
      break;
    case INPUT_VECTORISED:
      res[i] = dmpb_(x[i], alpha[i], beta[i], c[i]);
      break;
    }
  }

  if(log_p)
    res = log(res);

  return res;
}



//[[Rcpp::export]]
NumericVector cpp_pmpb(NumericVector q, NumericVector alpha, NumericVector beta, NumericVector c, const bool& lower_tail, const bool& log_p) {
  int n = q.size();
  NumericVector res(n);

  int max_q = max(q);
  if(1 == alpha.size() && 1 == beta.size() && 1 == c.size()) {
    // single parameters
    for(int i = 0; i < n; i++) {
      res[i] = pmpb_(q[i], alpha[0], beta[0], c[0]);
    }
  } else if(n == alpha.size() && n == beta.size() && n == c.size()) {
    for(int i = 0; i < n; i++) {
      res[i] = pmpb_(q[i], alpha[i], beta[i], c[i]);
    }
  } else {
    warning("Dimensions do not match");
  }

  if(!lower_tail)
    res = 1.0 - res;

  if(log_p)
    res = log(res);

  return res;
}


// [[Rcpp::export]]
NumericVector cpp_rmpb(double n, NumericVector alpha, NumericVector beta, NumericVector c) {
  if(isInadmissible(n)) {
    stop("Error in rmpb : invalid arguments");
  }

  NumericVector res((int)n);

  if(1 == alpha.size() && 1 == beta.size() && 1 == c.size()) {
    // single parameters
    if(isInadmissible(alpha[0]) || isInadmissible(beta[0]) || isInadmissible(c[0])) {
      NumericVector na_res((int)n, NA_REAL);
      warning("NAs produced");
      return na_res;
    }
    NumericVector poissonParameter = rbeta(n, alpha[0], beta[0]) * c[0];
    for(int i = 0; i < n; i++) {
      NumericVector t = rpois(1, poissonParameter[i]);
      res[i] = t[0];
    }
  } else if(n == alpha.size() && n == beta.size() && n == c.size()) {
    // vectorised parameters
    for(int i = 0; i < n; i++) {
      if(isInadmissible(alpha[i]) || isInadmissible(beta[i]) || isInadmissible(c[i])) {
        res[i] = NA_REAL;
        warning("NAs produced");
      } else {
        NumericVector poissonParameter = rbeta(1, alpha[i], beta[i]) * c[i];
        NumericVector t = rpois(1, poissonParameter[0]);
        res[i] = t[0];
      }
    }
  } else {
    warning("Dimensions do not match");
  }
  return res;
}


// [[Rcpp::export]]
NumericVector cpp_qmpb(NumericVector p, NumericVector alpha, NumericVector beta, NumericVector c, const bool& lower_tail, const bool& log_p) {
    int n = p.size();
    NumericVector res(n);

    if(log_p)
      p = exp(p);

    if(lower_tail)
      p = 1.0 - p;

    if (1 == alpha.size() && 1 == beta.size() && 1 == c.size()) {
      // single parameters
      if(isInadmissible(alpha[0]) || isInadmissible(beta[0]) || isInadmissible(c[0])) {
        NumericVector na_res(n, NA_REAL);
        return na_res;
      }

      double* p_distr = pmpb_(alpha[0], beta[0], c[0]);
      for(int i = 0; i < n; i++) {
        if(p[i] == 1.0) {
          res[i] = R_PosInf;
        } else {
          res[i] = qmpb_(p[i], p_distr);
          res[i] = res[i] == Q_LIMIT ? R_PosInf : res[i];
        }
      }
    } else if (n == alpha.size() && n == beta.size() && n == c.size()) {
      // vectorised parameters
      for(int i = 0; i < n; i++) {
        if(p[i] == 1.0) {
          res[i] = R_PosInf;
        } else {
          res[i] = qmpb_(p[i], alpha[i], beta[i], c[i]);
          res[i] = res[i] == Q_LIMIT ? R_PosInf : res[i];
        }
      }
    } else {
      warning("Dimensions do not match");
    }
    return res;
}
