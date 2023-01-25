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
double multinom_LL(double p, double sum_exp_p, bool z_ij, double k, double n, double mean, double precision) {
  double mm = 0;
  if(z_ij) {
    mm = k*p - n*log(exp(p) + sum_exp_p);
  }
  mm -= 0.5*precision*(p - mean)*(p - mean);
  return mm;
}


double cond_mv_mean(NumericVector x, NumericVector mean, NumericMatrix Q, int i) {
  double mm = mean[i];
  double Qii = Q(i, i);
  for(int j = 0; j < x.size(); j++) {
    if(i != j) {
      mm -= Q(i, j)*(x[j] - mean[j]) / Qii;
    }
  }
  return mm;
}

double cond_mv_mean_multinom(NumericVector x, IntegerVector refidx, NumericVector mean, NumericMatrix Q, int ij) {
  int i = refidx[ij];
  double mm = mean[i];
  double Qii = Q(i, i);
  for(int jj = 0; jj < x.size(); jj++) {
    int j = refidx[jj];
    if(ij != jj &&  j >= 0) {
      mm -= Q(i, j)*(x[jj] - mean[j]) / Qii;
    }
  }
  return mm;
}

