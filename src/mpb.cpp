// [[Rcpp::depends(RcppGSL)]]

#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include <Rcpp.h>
#include <cmath>
#include "mpb.h"
using namespace Rcpp;

// kummer series using GSL
double kummer_(double x, double a, double b, int lnchf) {
  double res = gsl_sf_hyperg_1F1(a, b, x);
  if(1 == lnchf) {
    return log(res);
  } else {
    return res;
  }
}

double dmpb_(double x, double alpha, double beta, double c) {
  // double cre, cim;
  // double are = alpha+x, aim = 0.0, bre = alpha+beta+x, bim = 0.0;
  // int n = 1, ip = 0, lnchf = 1;
  // double zre = -c;
  // double zim = 0;
  // chfm_(&zre, &zim, &are, &aim, &bre, &bim, &cre, &cim, &n, &lnchf, &ip);

  // using gsl
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
    // return denom;
  }
}

double pmpb_(double x, double alpha, double beta, double c) {
  double res = 0;
  for(int i = 0; i <= x; i++) {
    res += dmpb_(i, alpha, beta, c);
  }
  return res;
}

//'@rdname mpb2
//'@export
// [[Rcpp::export]]
NumericVector kummer_gsl(NumericVector x, double a, double b, int lnchf = 0) {
    int n = x.size();
    NumericVector res(n);
    for(int i = 0; i < n; i++) {
        res[i] = kummer_(x[i], a, b, lnchf);
    }
    return res;
}

//'@rdname mpb2
//'@export
//'@examples
//'  X <- dmpb(x=0:200, alpha=5, beta=3, c=20)
//'  plot(0:200, X, type='l')
//'  Y <- dmpb(0:10, seq(10.0,11.0,by=0.1), seq(30.0,31.0,by=0.1), seq(10.2,11.2,by=0.1))
// [[Rcpp::export]]
NumericVector dmpb(IntegerVector x, NumericVector alpha, NumericVector beta, NumericVector c) {
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
  return res;
}


//'@rdname mpb2
//'@export
//[[Rcpp::export]]
NumericVector pmpb(IntegerVector q, NumericVector alpha, NumericVector beta, NumericVector c) {
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
  return res;
}


//'@rdname mpb2
//'@export
//'@examples
//'  RV <- rmpb(n = 1000, alpha=5, beta= 3, c=20)
//'  plot(0 : 200, X, type="l")
//'  lines(density(RV), col="red")
//'  R2 <- rmpb(11, seq(10.0,11.0,by=0.1), seq(30.0,31.0,by=0.1), seq(10.2,11.2,by=0.1))
// [[Rcpp::export]]
IntegerVector rmpb(int n, NumericVector alpha, NumericVector beta, NumericVector c) {
    IntegerVector res(n);

    if(1 == alpha.size() && 1 == beta.size() && 1 == c.size()) {
      // single parameters
      NumericVector poissonParameter = rbeta(n, alpha[0], beta[0]) * c[0];
      for(int i = 0; i < n; i++) {
        NumericVector t = rpois(1, poissonParameter[i]);
        res[i] = t[0];
      }
    } else if(n == alpha.size() && n == beta.size() && n == c.size()) {
      // vectorised parameters
      for(int i = 0; i < n; i++) {
        NumericVector poissonParameter = rbeta(1, alpha[i], beta[i]) * c[i];
        NumericVector t = rpois(1, poissonParameter[0]);
        res[i] = t[0];
      }
    } else {
      warning("Dimensions do not match");
    }
    return res;
}
