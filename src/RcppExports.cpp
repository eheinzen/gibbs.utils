// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// impute_conj_mvnorm_mu_cpp
NumericMatrix impute_conj_mvnorm_mu_cpp(NumericMatrix y, NumericMatrix mu, LogicalMatrix impute, NumericMatrix Q, NumericVector mu0, NumericVector tau0);
RcppExport SEXP _gibbs_utils_impute_conj_mvnorm_mu_cpp(SEXP ySEXP, SEXP muSEXP, SEXP imputeSEXP, SEXP QSEXP, SEXP mu0SEXP, SEXP tau0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type impute(imputeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau0(tau0SEXP);
    rcpp_result_gen = Rcpp::wrap(impute_conj_mvnorm_mu_cpp(y, mu, impute, Q, mu0, tau0));
    return rcpp_result_gen;
END_RCPP
}
// mh_binom
NumericVector mh_binom(NumericVector p, NumericVector proposal, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision);
RcppExport SEXP _gibbs_utils_mh_binom(SEXP pSEXP, SEXP proposalSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type precision(precisionSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_binom(p, proposal, k, n, mean, precision));
    return rcpp_result_gen;
END_RCPP
}
// mh_binom_mv
NumericVector mh_binom_mv(NumericMatrix p, NumericMatrix proposal, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q);
RcppExport SEXP _gibbs_utils_mh_binom_mv(SEXP pSEXP, SEXP proposalSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_binom_mv(p, proposal, k, n, mean, Q));
    return rcpp_result_gen;
END_RCPP
}
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
double one_binom_slice_mv(NumericVector p, double k, double n, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_binom_slice_mv(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP iSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
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
NumericMatrix slice_sample_binom_mv(NumericMatrix p, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_binom_mv(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_binom_mv(p, k, n, mean, Q, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// one_multinom_slice_mv
double one_multinom_slice_mv(NumericVector p_j, LogicalVector z_j, double k, double n, NumericVector p_i, NumericVector mean, NumericMatrix Q, int i, int j, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_multinom_slice_mv(SEXP p_jSEXP, SEXP z_jSEXP, SEXP kSEXP, SEXP nSEXP, SEXP p_iSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP iSEXP, SEXP jSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p_j(p_jSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type z_j(z_jSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_i(p_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_multinom_slice_mv(p_j, z_j, k, n, p_i, mean, Q, i, j, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_multinom_mv
NumericMatrix slice_sample_multinom_mv(List p_j, LogicalMatrix z, NumericMatrix k, NumericMatrix n, NumericMatrix p_i, NumericMatrix mean, NumericMatrix Q, int j, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_multinom_mv(SEXP p_jSEXP, SEXP zSEXP, SEXP kSEXP, SEXP nSEXP, SEXP p_iSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP jSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type p_j(p_jSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_i(p_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_multinom_mv(p_j, z, k, n, p_i, mean, Q, j, w, nexpand, ncontract));
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
double one_pois_slice_mv(NumericVector L, double k, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_one_pois_slice_mv(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP iSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
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
NumericMatrix slice_sample_pois_mv(NumericMatrix L, NumericMatrix k, NumericMatrix mean, NumericMatrix Q, double w, int nexpand, int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_pois_mv(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_pois_mv(L, k, mean, Q, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gibbs_utils_impute_conj_mvnorm_mu_cpp", (DL_FUNC) &_gibbs_utils_impute_conj_mvnorm_mu_cpp, 6},
    {"_gibbs_utils_mh_binom", (DL_FUNC) &_gibbs_utils_mh_binom, 6},
    {"_gibbs_utils_mh_binom_mv", (DL_FUNC) &_gibbs_utils_mh_binom_mv, 6},
    {"_gibbs_utils_one_binom_slice", (DL_FUNC) &_gibbs_utils_one_binom_slice, 8},
    {"_gibbs_utils_slice_sample_binom", (DL_FUNC) &_gibbs_utils_slice_sample_binom, 8},
    {"_gibbs_utils_one_binom_slice_mv", (DL_FUNC) &_gibbs_utils_one_binom_slice_mv, 9},
    {"_gibbs_utils_slice_sample_binom_mv", (DL_FUNC) &_gibbs_utils_slice_sample_binom_mv, 8},
    {"_gibbs_utils_one_multinom_slice_mv", (DL_FUNC) &_gibbs_utils_one_multinom_slice_mv, 12},
    {"_gibbs_utils_slice_sample_multinom_mv", (DL_FUNC) &_gibbs_utils_slice_sample_multinom_mv, 11},
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
