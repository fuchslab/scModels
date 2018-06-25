#include <Rcpp.h>
#include <array>
#include "shared.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' Gillespie Algorithm for the poisson-beta
//'
//' Gillespie Algorithm to simulate from basic kinetic model of gene activation and mRNA transcription
//' @param n number of simulations
//' @param r_act polymerase binding rate / DNA activation rate
//' @param r_deact polymerase unbinding rate / DNA deactivation rate
//' @param r_on transcription rate
//' @param r_degr mRNA degradation rate
//' @name gmRNA
//' @rdname gmRNA
//' @export
// [[Rcpp::export]]
NumericVector gmRNA_switch(double n, double r_act, double r_deact, double r_on, double r_degr) {
  if(!isInteger(n)) {
    return NumericVector(0);
  }

  NumericVector res((int)n);
  double t0 = 0, tmax = 200, tx;
  int k;
  std::array<double, 3> x0{ {1, 0, 0} };
  std::array<double, 3>x;
  double r_act1, r_act2, r_act3, r_act4, r_actx;
  double tau, tau_stern, u;

  for(int i = 0; i < n; i++) {
    tx = t0;
    x = x0;
    while(tx < tmax) {
      // step 1
      r_act1 = r_act * x[0];
      r_act2 = r_deact * x[1];
      r_act3 = r_on * x[1];
      r_act4 = r_degr * x[2];
      r_actx = r_act1 + r_act2 + r_act3 + r_act4;

      // step 2
      NumericVector tau_vec = rexp(1, r_actx);
      tau = tau_vec[0];
      tau_stern = min(NumericVector::create(tau, tmax - tx));

      // step 3
      NumericVector u_vec = runif(1);
      u = u_vec[0];
      if(u <= r_act1/r_actx)
        k = 1;
      else if(u <= (r_act1+r_act2)/r_actx)
        k = 2;
      else if(u <= (r_act1+r_act2+r_act3)/r_actx)
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
