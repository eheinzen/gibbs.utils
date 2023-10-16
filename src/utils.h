
#ifndef GIBBS_UTILS_UTILS_H
#define GIBBS_UTILS_UTILS_H

#include <Rcpp.h>
double binom_LL(double p, double k, double n, double mean, double precision);
double pois_LL(double L, double k, double mean, double precision, double trunc_at, bool lower);
double multinom_LL(double p, double sum_exp_p, bool z_ij, double k, double n, double mean, double precision);


double cond_mv_mean(Rcpp::NumericVector x, Rcpp::NumericVector mean, Rcpp::NumericMatrix Q, int i);
double cond_mv_mean_multinom(Rcpp::NumericVector x, Rcpp::IntegerVector refidx, Rcpp::NumericVector mean, Rcpp::NumericMatrix Q, int ij);

bool accept_reject(double ratio);

#endif
