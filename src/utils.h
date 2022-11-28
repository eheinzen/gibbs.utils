
#ifndef GIBBS_UTILS_UTILS_H
#define GIBBS_UTILS_UTILS_H

#include <Rcpp.h>
double binom_LL(double p, double k, double n, double mean, double precision);
double pois_LL(double L, double k, double mean, double precision, double trunc_at, bool lower);
double multinom_LL(Rcpp::NumericVector p_ij, Rcpp::LogicalVector z_ij, Rcpp::IntegerVector which_i,
                   double k, double n, double mean, double precision, int ij);


Rcpp::NumericVector replace_it(Rcpp::NumericVector x, int i, double value);
double cond_mv_mean(Rcpp::NumericVector x, Rcpp::NumericVector mean, Rcpp::NumericMatrix Q, int i);
bool accept_reject(double ratio);

#endif
