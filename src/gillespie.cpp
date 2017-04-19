#include <Rcpp.h>
#include <array>
#include "shared.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' Gillespie Algorithm to simulate from basic kinetic model of gene activation and mRNA transcription
//' @param n number of simulations
//' @param lambda polymerase binding rate / DNA activation rate
//' @param gamma polymerase unbinding rate / DNA deactivation rate
//' @param r transcription rate
//' @param mu mRNA degradation rate
//' @name gmRNA
//' @rdname gmRNA
//' @export
// [[Rcpp::export]]
NumericVector gmRNA(double n, double lambda, double gamma, double r, double mu) {
  if(!isInteger(n)) {
    return NumericVector(0);
  }

  NumericVector res((int)n);
  double t0 = 0, tmax = 200, tx;
  int k;
  std::array<double, 3> x0{ {1, 0, 0} };
  std::array<double, 3>x;
  double lambda1, lambda2, lambda3, lambda4, lambdax;
  double tau, tau_stern, u;

  for(int i = 0; i < n; i++) {
    tx = t0;
    x = x0;
    while(tx < tmax) {
      // step 1
      lambda1 = lambda * x[0];
      lambda2 = gamma * x[1];
      lambda3 = r * x[1];
      lambda4 = mu * x[2];
      lambdax = lambda1 + lambda2 + lambda3 + lambda4;

      // step 2
      NumericVector tau_vec = rexp(1, lambdax);
      tau = tau_vec[0];
      tau_stern = min(NumericVector::create(tau, tmax - tx));

      // step 3
      NumericVector u_vec = runif(1);
      u = u_vec[0];
      if(u <= lambda1/lambdax)
        k = 1;
      else if(u <= (lambda1+lambda2)/lambdax)
        k = 2;
      else if(u <= (lambda1+lambda2+lambda3)/lambdax)
        k = 3;
      else
        k = 4;

      // step 4
      switch(k) {
      case 1:
        x[0]--;
        x[1]++;
        break;
      case 2:
        x[0]++;
        x[1]--;
        break;
      case 3:
        x[2]++;
        break;
      case 4:
        x[2]--;
        break;
      }

      // step 5
      tx = tx + tau_stern;
    }
    res[i] = x[2];
  }
  return res;
}
