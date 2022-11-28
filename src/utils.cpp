#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool accept_reject(double ratio) {
  return log(R::runif(0.0, 1.0)) <= ratio;
}


double binom_LL(double p, double k, double n, double mean, double precision) {
  return k*p - n*log(1.0 + exp(p)) - 0.5*precision*(p - mean)*(p - mean);
}
double pois_LL(double L, double k, double mean, double precision, double trunc_at, bool lower) {
  double eL = exp(L);
  double out = k*L - eL - 0.5*precision*(L - mean)*(L - mean);
  if(trunc_at >= 0.0) {
    out -= R::ppois(trunc_at, eL, lower, 1);
  }
  return out;
}
double multinom_LL(NumericVector p_ij, LogicalVector z_ij, IntegerVector which_i, double k, double n, double mean, double precision, int ij) {
  double mm = 0;
  if(z_ij[ij]) {
    // we still get the right log-likelihood if we include this when z_ij[ij] is false,
    // but the exp() costs time
    double nn = 0;
    int i = which_i[ij];
    for(int jj = 0; jj < p_ij.size(); jj++) {
      if(z_ij[jj] && which_i[jj] == i) {
        nn += exp(p_ij[jj]);
      }
    }
    mm = k*p_ij[ij] - n*log(nn);
  }
  mm -= 0.5*precision*(p_ij[ij] - mean)*(p_ij[ij] - mean);
  return mm;
}


NumericVector replace_it(NumericVector x, int i, double value) {
  NumericVector out = clone(x);
  out[i] = value;
  return out;
}

double cond_mv_mean(NumericVector x, NumericVector mean, NumericMatrix Q, int i) {
  double mm = mean[i];
  for(int j = 0; j < x.size(); j++) {
    if(i != j) {
      mm -= Q(i, j)*(x[j] - mean[j]) / Q(i, i);
    }
  }
  return mm;
}
