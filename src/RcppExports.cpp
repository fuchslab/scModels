// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_gmRNA_basic
NumericVector cpp_gmRNA_basic(double n, double r_on, double r_degr);
RcppExport SEXP _scModels_cpp_gmRNA_basic(SEXP nSEXP, SEXP r_onSEXP, SEXP r_degrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r_on(r_onSEXP);
    Rcpp::traits::input_parameter< double >::type r_degr(r_degrSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gmRNA_basic(n, r_on, r_degr));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gmRNA_switch
NumericVector cpp_gmRNA_switch(double n, double r_act, double r_deact, double r_on, double r_degr);
RcppExport SEXP _scModels_cpp_gmRNA_switch(SEXP nSEXP, SEXP r_actSEXP, SEXP r_deactSEXP, SEXP r_onSEXP, SEXP r_degrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r_act(r_actSEXP);
    Rcpp::traits::input_parameter< double >::type r_deact(r_deactSEXP);
    Rcpp::traits::input_parameter< double >::type r_on(r_onSEXP);
    Rcpp::traits::input_parameter< double >::type r_degr(r_degrSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gmRNA_switch(n, r_act, r_deact, r_on, r_degr));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gmRNA_burst
NumericVector cpp_gmRNA_burst(double n, double r_burst, double s_burst, double r_degr);
RcppExport SEXP _scModels_cpp_gmRNA_burst(SEXP nSEXP, SEXP r_burstSEXP, SEXP s_burstSEXP, SEXP r_degrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r_burst(r_burstSEXP);
    Rcpp::traits::input_parameter< double >::type s_burst(s_burstSEXP);
    Rcpp::traits::input_parameter< double >::type r_degr(r_degrSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gmRNA_burst(n, r_burst, s_burst, r_degr));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gmRNA_basic_burst
NumericVector cpp_gmRNA_basic_burst(double n, double r_on, double r_burst, double s_burst, double r_degr);
RcppExport SEXP _scModels_cpp_gmRNA_basic_burst(SEXP nSEXP, SEXP r_onSEXP, SEXP r_burstSEXP, SEXP s_burstSEXP, SEXP r_degrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r_on(r_onSEXP);
    Rcpp::traits::input_parameter< double >::type r_burst(r_burstSEXP);
    Rcpp::traits::input_parameter< double >::type s_burst(s_burstSEXP);
    Rcpp::traits::input_parameter< double >::type r_degr(r_degrSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gmRNA_basic_burst(n, r_on, r_burst, s_burst, r_degr));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rInvGaus
NumericVector cpp_rInvGaus(double n, double mu, double lambda);
RcppExport SEXP _scModels_cpp_rInvGaus(SEXP nSEXP, SEXP muSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rInvGaus(n, mu, lambda));
    return rcpp_result_gen;
END_RCPP
}
// chf_1F1
NumericVector chf_1F1(NumericVector x, NumericVector a, NumericVector b);
RcppExport SEXP _scModels_chf_1F1(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(chf_1F1(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// cpp_dpb
NumericVector cpp_dpb(NumericVector& x, NumericVector& alpha, NumericVector& beta, NumericVector& c, const bool& log_p);
RcppExport SEXP _scModels_cpp_dpb(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dpb(x, alpha, beta, c, log_p));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ppb
NumericVector cpp_ppb(NumericVector& q, NumericVector& alpha, NumericVector& beta, NumericVector& c, const bool& lower_tail, const bool& log_p);
RcppExport SEXP _scModels_cpp_ppb(SEXP qSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ppb(q, alpha, beta, c, lower_tail, log_p));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rpb
NumericVector cpp_rpb(const int& n, NumericVector& alpha, NumericVector& beta, NumericVector& c);
RcppExport SEXP _scModels_cpp_rpb(SEXP nSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rpb(n, alpha, beta, c));
    return rcpp_result_gen;
END_RCPP
}
// cpp_qpb
NumericVector cpp_qpb(NumericVector& p, NumericVector& alpha, NumericVector& beta, NumericVector& c, const bool& lower_tail, const bool& log_p);
RcppExport SEXP _scModels_cpp_qpb(SEXP pSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_qpb(p, alpha, beta, c, lower_tail, log_p));
    return rcpp_result_gen;
END_RCPP
}
