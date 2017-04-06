// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// chf_1F1_gsl
NumericVector chf_1F1_gsl(NumericVector x, NumericVector a, NumericVector b, const bool& log_v);
RcppExport SEXP mpb2_chf_1F1_gsl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP log_vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_v(log_vSEXP);
    rcpp_result_gen = Rcpp::wrap(chf_1F1_gsl(x, a, b, log_v));
    return rcpp_result_gen;
END_RCPP
}
// cpp_dmpb
NumericVector cpp_dmpb(NumericVector x, NumericVector alpha, NumericVector beta, NumericVector c, const bool& log_p);
RcppExport SEXP mpb2_cpp_dmpb(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dmpb(x, alpha, beta, c, log_p));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pmpb
NumericVector cpp_pmpb(NumericVector q, NumericVector alpha, NumericVector beta, NumericVector c, const bool& lower_tail, const bool& log_p);
RcppExport SEXP mpb2_cpp_pmpb(SEXP qSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pmpb(q, alpha, beta, c, lower_tail, log_p));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rmpb
NumericVector cpp_rmpb(double n, NumericVector alpha, NumericVector beta, NumericVector c);
RcppExport SEXP mpb2_cpp_rmpb(SEXP nSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rmpb(n, alpha, beta, c));
    return rcpp_result_gen;
END_RCPP
}
// cpp_qmpb
NumericVector cpp_qmpb(NumericVector p, NumericVector alpha, NumericVector beta, NumericVector c, const bool& lower_tail, const bool& log_p);
RcppExport SEXP mpb2_cpp_qmpb(SEXP pSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP lower_tailSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< const bool& >::type lower_tail(lower_tailSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_qmpb(p, alpha, beta, c, lower_tail, log_p));
    return rcpp_result_gen;
END_RCPP
}
