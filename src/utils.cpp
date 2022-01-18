#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double binom_LL(double p, double k, double n, double mean, double precision) {
  return k*p - n*log(1.0 + exp(p)) - 0.5*precision*(p - mean)*(p - mean);
}
double multinom_LL_mv(NumericVector p_j, LogicalVector z_j, double k, double n, NumericVector p_i, NumericVector mean, NumericMatrix Q, int i, int j) {
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
  for(int ii = 0; ii < p_i.size(); ii++) {
    mm -= 0.5*(1.0 + (i != ii))*(p_i[i] - mean[i])*Q(i, ii)*(p_i[ii] - mean[ii]);
  }
  return mm;
}



double pois_LL(double L, double k, double mean, double precision) {
  return k*L - exp(L) - 0.5*precision*(L - mean)*(L - mean);
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
