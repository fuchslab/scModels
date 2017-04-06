#include <Rcpp.h>
#include "shared.h"

bool validKummerParameters(double a, double b, bool warn) {
  if(a < 0) {
    if(warn) {
      char msg[55];
      std::sprintf(msg, "Parameter a is less than zero: %f", a);
      Rcpp::warning(msg);
    }
    return false;
  }
  if(b < a) {
    if(warn) {
      char msg[55];
      std::sprintf(msg, "Wrong parameters: b cannot be less than a: %f < %f", b, a);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::sprintf(msg, "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

bool validProbability(double p) {
  if (p >= 0.0 && p <= 1.0) {
    return true;
  } else {
    char msg[55];
    std::sprintf(msg, "invalid probability: %f", p);
    Rcpp::warning(msg);
    return false;
  }
}

bool isInadmissible(double x, bool warn) {
  if(Rcpp::NumericVector::is_na(x) || Rcpp::traits::is_nan<REALSXP>(x)) {
    if(warn)
      Rcpp::warning("NA/NaNs given in input");
    return true;
  } else {
    return false;
  }
}
