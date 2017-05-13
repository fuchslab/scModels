// [[Rcpp::depends(RcppGSL)]]

#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_hyperg.h>
#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include "shared.h"
using namespace Rcpp;

// kummer series using GSL
double kummer_(double x, double a, double b, bool log_v) {
  if(!validKummerParameters(a, b)) {
    return R_NaN;
  }

  gsl_set_error_handler_off();
  gsl_sf_result gsl_res;
  double log_chf;
  int status = gsl_sf_hyperg_1F1_e(a, b, x, &gsl_res);
  if( status ) {
    if( status == GSL_EUNDRFLW || status == GSL_EOVRFLW ){
      int status_transform = gsl_sf_hyperg_1F1_e(b-a, b, -x, &gsl_res);
      if(status_transform) {
        warning("Kummer transformation failed!");
        return R_NaN;
      } else {
        log_chf = x + log(gsl_res.val);
        warning("using transformation");
      }
    } else {
      reportGslError(status);
      log_chf = 1;
    }
  } else {
    log_chf = log(gsl_res.val);
  }

  if(log_v) {
    return log_chf;
  } else {
    return std::exp(log_chf);
  }
}

// density function
double dmpb_(double x, double alpha, double beta, double c, bool& throw_warning) {
  if( isInadmissible(x) || isInadmissible(alpha) || isInadmissible(beta) || isInadmissible(c) )
    return x+alpha+beta+c;

  if( !isInteger(x) || x < 0  || traits::is_infinite<REALSXP>(x) )
    return 0;

  if(!validMpbParameters(alpha, beta, c)) {
    throw_warning = true;
    return R_NaN;
  }

  double cre = kummer_(-c, alpha+x, beta+alpha+x, true);
  if(isInadmissible(cre))
    return R_NaN;

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
double pmpb_(double x, double alpha, double beta, double c, bool& throw_warning) {
  if( isInadmissible(x) || isInadmissible(alpha) || isInadmissible(beta) || isInadmissible(c) )
    return x+alpha+beta+c;

  if(!validMpbParameters(alpha, beta, c)) {
    throw_warning = true;
    return R_NaN;
  }

  if( !isInteger(x) )
    return 0;
  if(traits::is_infinite<REALSXP>(x))
    return 1;
  double res = 0;
  for(int i = 0; i <= x; i++) {
    res += dmpb_(i, alpha, beta, c, throw_warning);
  }
  return res;
}

// distribution function array
double* pmpb_(double alpha, double beta, double c) {
  double *res = (double *)std::malloc(Q_LIMIT * sizeof(double));
  bool throw_warning = false;
  res[0] = dmpb_(0, alpha, beta, c, throw_warning);
  for(int i = 1; i < Q_LIMIT; i++) {
    res[i] = res[i-1] + dmpb_(i, alpha, beta, c, throw_warning);
  }
  return res;
}

// quantiles for single parameters
double qmpb_(double p, double *p_distr) {
  if(isInadmissible(p))
    return NA_REAL;
  if(!validProbability(p) || isInadmissible(p_distr[0])){
    warning("NaNs produced");
    return R_NaN;
  }

  if(p == 0.0)
    return 0.0;
  if(p == 1.0 || p > p_distr[Q_LIMIT-1])
    return R_PosInf;

  for(int i = 1; i < Q_LIMIT; i++) {
    if(p > p_distr[i-1] && p < p_distr[i]) {
      return i;
    }
  }

  return R_PosInf;
}

// quantiles for vectorised parameters
double qmpb_(double p, double alpha, double beta, double c) {
  if(isInadmissible(p) || isInadmissible(alpha) || isInadmissible(beta) || isInadmissible(c))
    return NA_REAL;
  if(!validProbability(p)){
    warning("NaNs produced");
    return R_NaN;
  }

  if(p == 0.0)
    return 0.0;

  double *p_distr = pmpb_(alpha, beta, c);

  if(p == 1.0 || p > p_distr[Q_LIMIT-1])
    return R_PosInf;

  for(int i = 1; i < Q_LIMIT; i++) {
    if(p > p_distr[i-1] && p < p_distr[i]) {
      return i;
    }
  }

  return R_PosInf;
}

// random number generator
double rmpb_(double alpha, double beta, double c, bool& throw_warning) {
  if(isInadmissible(alpha) || isInadmissible(beta) || isInadmissible(c)) {
    throw_warning = true;
    return NA_REAL;
  }

  if(!validMpbParameters(alpha, beta, c)) {
    throw_warning = true;
    return R_NaN;
  }

  NumericVector poissonParameter = rbeta(1, alpha, beta) * c;
  NumericVector t = rpois(1, poissonParameter[0]);

  return t[0];
}

//' Kummer's (confluent hypergeometric) function
//'
//' Kummer's (confluent hypergeometric) function of the first kind
//' for numeric (non-complex) values and input parameters
//' @param x numeric value or vector
//' @param a,b numeric parameters of the Kummer function
//' @param log_v logical; if TRUE, the log of the value is returned
//' @name chf_1F1_gsl
//' @rdname chf_1F1_gsl
//' @importFrom RcppGSL CFlags LdFlags
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
NumericVector cpp_dmpb(NumericVector& x, NumericVector& alpha, NumericVector& beta, NumericVector& c, const bool& log_p = false) {
  if(std::min({x.length(), alpha.length(), beta.length(), c.length()}) < 1) {
    return NumericVector(0);
  }

  int n = std::max({x.length(), alpha.length(), beta.length(), c.length()});
  NumericVector p(n);
  bool throw_warning = false;

  for(int i = 0; i < n; i++) {
    p[i] = dmpb_(GETV(x, i), GETV(alpha, i), GETV(beta, i), GETV(c, i), throw_warning);
  }

  if(log_p)
    p = log(p);

  if(throw_warning)
    warning("NaNs produced");

  return p;
}



//[[Rcpp::export]]
NumericVector cpp_pmpb(NumericVector& q, NumericVector& alpha, NumericVector& beta, NumericVector& c, const bool& lower_tail, const bool& log_p) {
  if(std::min({ q.length(), alpha.length(), beta.length(), c.length() }) < 1) {
    return NumericVector(0);
  }

  int n = std::max({ q.length(), alpha.length(), beta.length(), c.length() });
  NumericVector p(n);

  bool throw_warning = false;

  for(int i = 0; i < n; i++) {
    p[i] = pmpb_(GETV(q, i), GETV(alpha, i), GETV(beta, i), GETV(c, i), throw_warning);
  }

  if(!lower_tail)
    p = 1.0 - p;

  if(log_p)
    p = log(p);

  if(throw_warning)
    warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rmpb(const int& n, NumericVector& alpha, NumericVector& beta, NumericVector& c) {
  if(std::min({ alpha.length(), beta.length(), c.length() }) < 1) {
    warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  bool throw_warning = false;

  for(int i = 0; i < n; i++) {
    x[i] = rmpb_(GETV(alpha, i), GETV(beta, i), GETV(c, i), throw_warning);
  }

  if(throw_warning)
    warning("NAs produced");

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_qmpb(NumericVector& p, NumericVector& alpha, NumericVector& beta, NumericVector& c, const bool& lower_tail, const bool& log_p) {
  if(std::min({ p.length(), alpha.length(), beta.length(), c.length() }) < 1) {
    return NumericVector(0);
  }

  int n = std::max({ p.length(), alpha.length(), beta.length(), c.length()});
  NumericVector res(n);

  if(log_p)
    p = exp(p);

  if(lower_tail)
    p = 1.0 - p;

  if (min(alpha) == max(alpha) && min(beta) == max(beta) && min(c) == max(c)) {
    // single parameters
    // optmized to compute cdf only once
    if(isInadmissible(alpha[0]) || isInadmissible(beta[0]) || isInadmissible(c[0])) {
      return NumericVector(n, NA_REAL);
    } else {
      double* p_distr = pmpb_(min(na_omit(alpha)), min(na_omit(beta)), min(na_omit(c)));
      for(int i = 0; i < n; i++) {
        res[i] = qmpb_(GETV(p, i), p_distr);
      }
    }
  } else {
    // vectorised parameters
    for(int i = 0; i < n; i++) {
      res[i] = qmpb_(GETV(p, i), GETV(alpha, i), GETV(beta, i), GETV(c, i));
    }
  }
  return res;
}
