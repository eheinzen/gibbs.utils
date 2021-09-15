#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double binom_LL(double p, double k, double n, double mean, double precision) {
  return k*p - n*log(1.0 + exp(p)) - 0.5*precision*(p - mean)*(p - mean);
}
double binom_LL_mv(NumericVector p, double k, double n, NumericVector mean, NumericMatrix Q, int i) {
  double mm = k*p[i] - n*log(1.0 + exp(p[i])) ;
  for(int j = 0; j < p.size(); j++) {
    mm -= 0.5*(1.0 + (i != j))*(p[i] - mean[i])*Q(i, j)*(p[j] - mean[j]);
  }
  return mm;
}


double multinom_LL_mv(NumericVector p_j, LogicalVector z_j, double k, double n, NumericVector p_i, NumericVector mean, NumericMatrix Q, int i, int j) {
  double nn = 0;
  for(int jj = 0; jj < p_j.size(); jj++) {
    if(z_j[jj]) {
      nn += exp(p_j[jj]);
    }
  }
  double mm = k*p_j[j] - n*log(nn);
  for(int ii = 0; ii < p_i.size(); ii++) {
    mm -= 0.5*(1.0 + (i != ii))*(p_i[i] - mean[i])*Q(i, ii)*(p_i[ii] - mean[ii]);
  }
  return mm;
}



double pois_LL(double L, double k, double mean, double precision) {
  return k*L - exp(L) - 0.5*precision*(L - mean)*(L - mean);
}
double pois_LL_mv(NumericVector L, double k, NumericVector mean, NumericMatrix Q, int i) {
  double mm = k*L[i] - exp(L[i]);
  for(int j = 0; j < L.size(); j++) {
    mm -= 0.5*(1.0 + (i != j))*(L[i] - mean[i])*Q(i, j)*(L[j] - mean[j]);
  }
  return mm;
}

NumericVector replace_it(NumericVector x, int i, double value) {
  NumericVector out = clone(x);
  out[i] = value;
  return out;
}
