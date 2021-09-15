
#ifndef GIBBS_UTILS_UTILS_H
#define GIBBS_UTILS_UTILS_H

#include <Rcpp.h>
double binom_LL(double p, double k, double n, double mean, double precision);
double binom_LL_mv(Rcpp::NumericVector p, double k, double n,
                   Rcpp::NumericVector mean, Rcpp::NumericMatrix Q, int i);
double pois_LL(double L, double k, double mean, double precision);
double pois_LL_mv(Rcpp::NumericVector L, double k, Rcpp::NumericVector mean,
                  Rcpp::NumericMatrix Q, int i);

Rcpp::NumericVector replace_it(Rcpp::NumericVector x, int i, double value);

#endif
