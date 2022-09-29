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
double multinom_LL(NumericVector p_j, LogicalVector z_j, double k, double n, double mean, double precision, int j) {
  double mm = 0;
  if(z_j[j]) {
    // we still get the right log-likelihood if we include this when z_j[j] is false,
    // but the exp() costs time
    double nn = 0;
    for(int jj = 0; jj < p_j.size(); jj++) {
      if(z_j[jj]) {
        nn += exp(p_j[jj]);
      }
    }
    mm = k*p_j[j] - n*log(nn);
  }
  mm -= 0.5*precision*(p_j[j] - mean)*(p_j[j] - mean);
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
