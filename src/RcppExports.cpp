// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// impute_conj_mvnorm_mu_cpp
NumericMatrix impute_conj_mvnorm_mu_cpp(const NumericMatrix y, const NumericMatrix mu, const LogicalMatrix impute, const NumericMatrix Q, const NumericVector mu0, const NumericVector tau0);
RcppExport SEXP _gibbs_utils_impute_conj_mvnorm_mu_cpp(SEXP ySEXP, SEXP muSEXP, SEXP imputeSEXP, SEXP QSEXP, SEXP mu0SEXP, SEXP tau0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type impute(imputeSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type tau0(tau0SEXP);
    rcpp_result_gen = Rcpp::wrap(impute_conj_mvnorm_mu_cpp(y, mu, impute, Q, mu0, tau0));
    return rcpp_result_gen;
END_RCPP
}
// mh_binom
NumericVector mh_binom(const bool qt, const NumericVector p, const NumericVector proposal, const NumericVector k, const NumericVector n, const NumericVector mean, const NumericVector precision, const int acceptance);
RcppExport SEXP _gibbs_utils_mh_binom(SEXP qtSEXP, SEXP pSEXP, SEXP proposalSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP acceptanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool >::type qt(qtSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const int >::type acceptance(acceptanceSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_binom(qt, p, proposal, k, n, mean, precision, acceptance));
    return rcpp_result_gen;
END_RCPP
}
// mh_binom_mv
NumericVector mh_binom_mv(const bool qt, const NumericMatrix p, const NumericMatrix proposal, const NumericMatrix k, const NumericMatrix n, const NumericMatrix mean, const NumericMatrix Q, const LogicalVector use_norm, const NumericMatrix norm, const int acceptance);
RcppExport SEXP _gibbs_utils_mh_binom_mv(SEXP qtSEXP, SEXP pSEXP, SEXP proposalSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP use_normSEXP, SEXP normSEXP, SEXP acceptanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool >::type qt(qtSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type use_norm(use_normSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type norm(normSEXP);
    Rcpp::traits::input_parameter< const int >::type acceptance(acceptanceSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_binom_mv(qt, p, proposal, k, n, mean, Q, use_norm, norm, acceptance));
    return rcpp_result_gen;
END_RCPP
}
// mh_multinom_mv
NumericVector mh_multinom_mv(const bool qt, const NumericMatrix p_ij, const NumericMatrix proposal, const LogicalMatrix z_ij, const IntegerVector which_i, const LogicalVector is_ref, const NumericMatrix k_ij, const NumericMatrix n_ij, const NumericMatrix mean, const NumericMatrix Q, const NumericVector Qdiag, const LogicalVector use_norm, const NumericMatrix norm, const bool diag, const int acceptance);
RcppExport SEXP _gibbs_utils_mh_multinom_mv(SEXP qtSEXP, SEXP p_ijSEXP, SEXP proposalSEXP, SEXP z_ijSEXP, SEXP which_iSEXP, SEXP is_refSEXP, SEXP k_ijSEXP, SEXP n_ijSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP QdiagSEXP, SEXP use_normSEXP, SEXP normSEXP, SEXP diagSEXP, SEXP acceptanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool >::type qt(qtSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type p_ij(p_ijSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type z_ij(z_ijSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type which_i(which_iSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type is_ref(is_refSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type k_ij(k_ijSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type n_ij(n_ijSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type Qdiag(QdiagSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type use_norm(use_normSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type norm(normSEXP);
    Rcpp::traits::input_parameter< const bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< const int >::type acceptance(acceptanceSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_multinom_mv(qt, p_ij, proposal, z_ij, which_i, is_ref, k_ij, n_ij, mean, Q, Qdiag, use_norm, norm, diag, acceptance));
    return rcpp_result_gen;
END_RCPP
}
// mh_pois
NumericVector mh_pois(const int method, const NumericVector L, const NumericVector proposal, const NumericVector k, const LogicalVector k_na, const NumericVector mean, const LogicalVector mean_inf, const NumericVector precision, const NumericVector trunc_at, const LogicalVector lower, const int acceptance);
RcppExport SEXP _gibbs_utils_mh_pois(SEXP methodSEXP, SEXP LSEXP, SEXP proposalSEXP, SEXP kSEXP, SEXP k_naSEXP, SEXP meanSEXP, SEXP mean_infSEXP, SEXP precisionSEXP, SEXP trunc_atSEXP, SEXP lowerSEXP, SEXP acceptanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type k_na(k_naSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type mean_inf(mean_infSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type trunc_at(trunc_atSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const int >::type acceptance(acceptanceSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_pois(method, L, proposal, k, k_na, mean, mean_inf, precision, trunc_at, lower, acceptance));
    return rcpp_result_gen;
END_RCPP
}
// mh_pois_mv
NumericVector mh_pois_mv(const int method, const NumericMatrix L, const NumericMatrix proposal, const NumericMatrix k, const LogicalMatrix k_na, const NumericMatrix mean, const NumericMatrix Q, const NumericMatrix trunc_at, const LogicalMatrix lower, const LogicalVector use_norm, const NumericMatrix norm, const int acceptance);
RcppExport SEXP _gibbs_utils_mh_pois_mv(SEXP methodSEXP, SEXP LSEXP, SEXP proposalSEXP, SEXP kSEXP, SEXP k_naSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP trunc_atSEXP, SEXP lowerSEXP, SEXP use_normSEXP, SEXP normSEXP, SEXP acceptanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type L(LSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type k_na(k_naSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type trunc_at(trunc_atSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type use_norm(use_normSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type norm(normSEXP);
    Rcpp::traits::input_parameter< const int >::type acceptance(acceptanceSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_pois_mv(method, L, proposal, k, k_na, mean, Q, trunc_at, lower, use_norm, norm, acceptance));
    return rcpp_result_gen;
END_RCPP
}
// one_binom_slice
double one_binom_slice(const double p, const double k, const double n, const double mean, const double precision, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_one_binom_slice(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_binom_slice(p, k, n, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_binom
NumericVector slice_sample_binom(const NumericVector p, const NumericVector k, const NumericVector n, const NumericVector mean, const NumericVector precision, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_binom(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_binom(p, k, n, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_binom_mv
NumericMatrix slice_sample_binom_mv(const NumericMatrix p, const NumericMatrix k, const NumericMatrix n, const NumericMatrix mean, const NumericMatrix Q, const LogicalVector use_norm, const NumericMatrix norm, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_binom_mv(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP use_normSEXP, SEXP normSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type use_norm(use_normSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type norm(normSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_binom_mv(p, k, n, mean, Q, use_norm, norm, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// one_multinom_slice
double one_multinom_slice(const double p_ij, const double sum_exp_p, const bool z_ij, const double k, const double n, const double mean, const double precision, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_one_multinom_slice(SEXP p_ijSEXP, SEXP sum_exp_pSEXP, SEXP z_ijSEXP, SEXP kSEXP, SEXP nSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type p_ij(p_ijSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_exp_p(sum_exp_pSEXP);
    Rcpp::traits::input_parameter< const bool >::type z_ij(z_ijSEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_multinom_slice(p_ij, sum_exp_p, z_ij, k, n, mean, precision, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_multinom_mv
NumericMatrix slice_sample_multinom_mv(const NumericMatrix p_ij, const LogicalMatrix z_ij, const IntegerVector which_i, const LogicalVector is_ref, const NumericMatrix k_ij, const NumericMatrix n_ij, const NumericMatrix mean, const NumericMatrix Q, const NumericVector Qdiag, const bool diag, const LogicalVector use_norm, const NumericMatrix norm, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_multinom_mv(SEXP p_ijSEXP, SEXP z_ijSEXP, SEXP which_iSEXP, SEXP is_refSEXP, SEXP k_ijSEXP, SEXP n_ijSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP QdiagSEXP, SEXP diagSEXP, SEXP use_normSEXP, SEXP normSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type p_ij(p_ijSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type z_ij(z_ijSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type which_i(which_iSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type is_ref(is_refSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type k_ij(k_ijSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type n_ij(n_ijSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type Qdiag(QdiagSEXP);
    Rcpp::traits::input_parameter< const bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type use_norm(use_normSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type norm(normSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_multinom_mv(p_ij, z_ij, which_i, is_ref, k_ij, n_ij, mean, Q, Qdiag, diag, use_norm, norm, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// one_pois_slice
double one_pois_slice(const double L, const double k, const double mean, const double precision, const double trunc_at, const bool lower, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_one_pois_slice(SEXP LSEXP, SEXP kSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP trunc_atSEXP, SEXP lowerSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type L(LSEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const double >::type trunc_at(trunc_atSEXP);
    Rcpp::traits::input_parameter< const bool >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(one_pois_slice(L, k, mean, precision, trunc_at, lower, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_pois
NumericVector slice_sample_pois(const NumericVector L, const NumericVector k, const LogicalVector k_na, const NumericVector mean, const LogicalVector mean_inf, const NumericVector precision, const NumericVector trunc_at, const LogicalVector lower, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_pois(SEXP LSEXP, SEXP kSEXP, SEXP k_naSEXP, SEXP meanSEXP, SEXP mean_infSEXP, SEXP precisionSEXP, SEXP trunc_atSEXP, SEXP lowerSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type k_na(k_naSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type mean_inf(mean_infSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type trunc_at(trunc_atSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_pois(L, k, k_na, mean, mean_inf, precision, trunc_at, lower, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// slice_sample_pois_mv
NumericMatrix slice_sample_pois_mv(const NumericMatrix L, const NumericMatrix k, const LogicalMatrix k_na, const NumericMatrix mean, const NumericMatrix Q, const NumericMatrix trunc_at, const LogicalMatrix lower, const LogicalVector use_norm, const NumericMatrix norm, const double w, const int nexpand, const int ncontract);
RcppExport SEXP _gibbs_utils_slice_sample_pois_mv(SEXP LSEXP, SEXP kSEXP, SEXP k_naSEXP, SEXP meanSEXP, SEXP QSEXP, SEXP trunc_atSEXP, SEXP lowerSEXP, SEXP use_normSEXP, SEXP normSEXP, SEXP wSEXP, SEXP nexpandSEXP, SEXP ncontractSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type L(LSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type k_na(k_naSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type trunc_at(trunc_atSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type use_norm(use_normSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type norm(normSEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nexpand(nexpandSEXP);
    Rcpp::traits::input_parameter< const int >::type ncontract(ncontractSEXP);
    rcpp_result_gen = Rcpp::wrap(slice_sample_pois_mv(L, k, k_na, mean, Q, trunc_at, lower, use_norm, norm, w, nexpand, ncontract));
    return rcpp_result_gen;
END_RCPP
}
// accept_reject
bool accept_reject(const double ratio);
RcppExport SEXP _gibbs_utils_accept_reject(SEXP ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type ratio(ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(accept_reject(ratio));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gibbs_utils_impute_conj_mvnorm_mu_cpp", (DL_FUNC) &_gibbs_utils_impute_conj_mvnorm_mu_cpp, 6},
    {"_gibbs_utils_mh_binom", (DL_FUNC) &_gibbs_utils_mh_binom, 8},
    {"_gibbs_utils_mh_binom_mv", (DL_FUNC) &_gibbs_utils_mh_binom_mv, 10},
    {"_gibbs_utils_mh_multinom_mv", (DL_FUNC) &_gibbs_utils_mh_multinom_mv, 15},
    {"_gibbs_utils_mh_pois", (DL_FUNC) &_gibbs_utils_mh_pois, 11},
    {"_gibbs_utils_mh_pois_mv", (DL_FUNC) &_gibbs_utils_mh_pois_mv, 12},
    {"_gibbs_utils_one_binom_slice", (DL_FUNC) &_gibbs_utils_one_binom_slice, 8},
    {"_gibbs_utils_slice_sample_binom", (DL_FUNC) &_gibbs_utils_slice_sample_binom, 8},
    {"_gibbs_utils_slice_sample_binom_mv", (DL_FUNC) &_gibbs_utils_slice_sample_binom_mv, 10},
    {"_gibbs_utils_one_multinom_slice", (DL_FUNC) &_gibbs_utils_one_multinom_slice, 10},
    {"_gibbs_utils_slice_sample_multinom_mv", (DL_FUNC) &_gibbs_utils_slice_sample_multinom_mv, 15},
    {"_gibbs_utils_one_pois_slice", (DL_FUNC) &_gibbs_utils_one_pois_slice, 9},
    {"_gibbs_utils_slice_sample_pois", (DL_FUNC) &_gibbs_utils_slice_sample_pois, 11},
    {"_gibbs_utils_slice_sample_pois_mv", (DL_FUNC) &_gibbs_utils_slice_sample_pois_mv, 12},
    {"_gibbs_utils_accept_reject", (DL_FUNC) &_gibbs_utils_accept_reject, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_gibbs_utils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
