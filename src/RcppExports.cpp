// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// one_binom_slice
double one_binom_slice(double p, double k, double n, double mean, double precision, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_binom_slice(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_binom_slice(p, k, n, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_binom
NumericVector slice_sample_binom(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_binom(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_binom(p, k, n, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// one_binom_slice_mv
double one_binom_slice_mv(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_binom_slice_mv(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP iSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_binom_slice_mv(p, k, n, mean, Q, i, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_binom_mv
NumericVector slice_sample_binom_mv(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericMatrix Q, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_binom_mv(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_binom_mv(p, k, n, mean, Q, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// one_pois_slice
double one_pois_slice(double L, double k, double mean, double precision, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_pois_slice(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_pois_slice(L, k, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_pois
NumericVector slice_sample_pois(NumericVector L, NumericVector k, NumericVector mean, NumericVector precision, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_pois(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_pois(L, k, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// one_pois_slice_mv
double one_pois_slice_mv(NumericVector L, NumericVector k, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_pois_slice_mv(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP iSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_pois_slice_mv(L, k, mean, Q, i, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_pois_mv
NumericVector slice_sample_pois_mv(NumericVector L, NumericVector k, NumericVector mean, NumericMatrix Q, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_pois_mv(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_pois_mv(L, k, mean, Q, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gibbs_utils_one_binom_slice", (DL_FUNC) &_gibbs_utils_one_binom_slice, 8},
    {"_gibbs_utils_slice_sample_binom", (DL_FUNC) &_gibbs_utils_slice_sample_binom, 8},
    {"_gibbs_utils_one_binom_slice_mv", (DL_FUNC) &_gibbs_utils_one_binom_slice_mv, 9},
    {"_gibbs_utils_slice_sample_binom_mv", (DL_FUNC) &_gibbs_utils_slice_sample_binom_mv, 8},
    {"_gibbs_utils_one_pois_slice", (DL_FUNC) &_gibbs_utils_one_pois_slice, 7},
    {"_gibbs_utils_slice_sample_pois", (DL_FUNC) &_gibbs_utils_slice_sample_pois, 7},
    {"_gibbs_utils_one_pois_slice_mv", (DL_FUNC) &_gibbs_utils_one_pois_slice_mv, 8},
    {"_gibbs_utils_slice_sample_pois_mv", (DL_FUNC) &_gibbs_utils_slice_sample_pois_mv, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_gibbs_utils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
