# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

impute_conj_mvnorm_mu_cpp <- function(y, mu, impute, Q, mu0, tau0) {
    .Call('_gibbs_utils_impute_conj_mvnorm_mu_cpp', PACKAGE = 'gibbs.utils', y, mu, impute, Q, mu0, tau0)
}

m_binom <- function(p, proposal, k, n, mean, precision) {
    .Call('_gibbs_utils_m_binom', PACKAGE = 'gibbs.utils', p, proposal, k, n, mean, precision)
}

m_binom_mv <- function(p, proposal, k, n, mean, Q, use_norm, norm) {
    .Call('_gibbs_utils_m_binom_mv', PACKAGE = 'gibbs.utils', p, proposal, k, n, mean, Q, use_norm, norm)
}

qt_binom_approx <- function(around, k, n, mean, tau) {
    .Call('_gibbs_utils_qt_binom_approx', PACKAGE = 'gibbs.utils', around, k, n, mean, tau)
}

qt_binom <- function(p, k, n, mean, precision) {
    .Call('_gibbs_utils_qt_binom', PACKAGE = 'gibbs.utils', p, k, n, mean, precision)
}

qt_binom_mv <- function(p, k, n, mean, Q, use_norm, norm) {
    .Call('_gibbs_utils_qt_binom_mv', PACKAGE = 'gibbs.utils', p, k, n, mean, Q, use_norm, norm)
}

m_pois <- function(L, proposal, k, k_na, mean, precision) {
    .Call('_gibbs_utils_m_pois', PACKAGE = 'gibbs.utils', L, proposal, k, k_na, mean, precision)
}

m_pois_mv <- function(L, proposal, k, k_na, mean, Q, use_norm, norm) {
    .Call('_gibbs_utils_m_pois_mv', PACKAGE = 'gibbs.utils', L, proposal, k, k_na, mean, Q, use_norm, norm)
}

qt_pois_approx <- function(around, k, mean, tau) {
    .Call('_gibbs_utils_qt_pois_approx', PACKAGE = 'gibbs.utils', around, k, mean, tau)
}

qt_pois <- function(L, k, k_na, mean, precision) {
    .Call('_gibbs_utils_qt_pois', PACKAGE = 'gibbs.utils', L, k, k_na, mean, precision)
}

qt_pois_mv <- function(L, k, k_na, mean, Q, use_norm, norm) {
    .Call('_gibbs_utils_qt_pois_mv', PACKAGE = 'gibbs.utils', L, k, k_na, mean, Q, use_norm, norm)
}

one_binom_slice <- function(p, k, n, mean, precision, w, nexpand, ncontract) {
    .Call('_gibbs_utils_one_binom_slice', PACKAGE = 'gibbs.utils', p, k, n, mean, precision, w, nexpand, ncontract)
}

slice_sample_binom <- function(p, k, n, mean, precision, w, nexpand, ncontract) {
    .Call('_gibbs_utils_slice_sample_binom', PACKAGE = 'gibbs.utils', p, k, n, mean, precision, w, nexpand, ncontract)
}

slice_sample_binom_mv <- function(p, k, n, mean, Q, use_norm, norm, w, nexpand, ncontract) {
    .Call('_gibbs_utils_slice_sample_binom_mv', PACKAGE = 'gibbs.utils', p, k, n, mean, Q, use_norm, norm, w, nexpand, ncontract)
}

one_multinom_slice <- function(p_j, z_j, k, n, mean, precision, j, w, nexpand, ncontract) {
    .Call('_gibbs_utils_one_multinom_slice', PACKAGE = 'gibbs.utils', p_j, z_j, k, n, mean, precision, j, w, nexpand, ncontract)
}

slice_sample_multinom_mv <- function(p_j, z, k, n, p_i, mean, Q, j, w, nexpand, ncontract) {
    .Call('_gibbs_utils_slice_sample_multinom_mv', PACKAGE = 'gibbs.utils', p_j, z, k, n, p_i, mean, Q, j, w, nexpand, ncontract)
}

one_pois_slice <- function(L, k, mean, precision, w, nexpand, ncontract) {
    .Call('_gibbs_utils_one_pois_slice', PACKAGE = 'gibbs.utils', L, k, mean, precision, w, nexpand, ncontract)
}

slice_sample_pois <- function(L, k, k_na, mean, precision, w, nexpand, ncontract) {
    .Call('_gibbs_utils_slice_sample_pois', PACKAGE = 'gibbs.utils', L, k, k_na, mean, precision, w, nexpand, ncontract)
}

slice_sample_pois_mv <- function(L, k, k_na, mean, Q, use_norm, norm, w, nexpand, ncontract) {
    .Call('_gibbs_utils_slice_sample_pois_mv', PACKAGE = 'gibbs.utils', L, k, k_na, mean, Q, use_norm, norm, w, nexpand, ncontract)
}

accept_reject <- function(ratio) {
    .Call('_gibbs_utils_accept_reject', PACKAGE = 'gibbs.utils', ratio)
}

