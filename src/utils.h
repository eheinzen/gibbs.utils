
#ifndef GIBBS_UTILS_UTILS_H
#define GIBBS_UTILS_UTILS_H

#include <Rcpp.h>
double binom_LL(double p, double k, double n, double mean, double precision);
double pois_LL(double L, double k, double mean, double precision);
double multinom_LL(Rcpp::NumericVector p_j, Rcpp::LogicalVector z_j, double k, double n,
                      double mean, double precision, int j);


Rcpp::NumericVector replace_it(Rcpp::NumericVector x, int i, double value);
double cond_mv_mean(Rcpp::NumericVector x, Rcpp::NumericVector mean, Rcpp::NumericMatrix Q, int i);
bool accept_reject(double ratio);

#endif
