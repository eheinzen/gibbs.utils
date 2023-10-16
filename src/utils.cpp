#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool accept_reject(const double ratio) {
  return log(R::runif(0.0, 1.0)) <= ratio;
}


double binom_LL(const double p, const double k, const double n,
                const double mean, const double precision) {
  return k*p - n*log(1.0 + exp(p)) - 0.5*precision*(p - mean)*(p - mean);
}
double pois_LL(const double L, const double k,
               const double mean, const double precision,
               const double trunc_at, const bool lower) {
  double eL = exp(L);
  double out = k*L - eL - 0.5*precision*(L - mean)*(L - mean);
  if(trunc_at >= 0.0) {
    out -= R::ppois(trunc_at, eL, lower, 1);
  }
  return out;
}
double multinom_LL(const double p, const double sum_exp_p, const bool z_ij,
                   const double k, const double n,
                   const double mean, const double precision) {
  double mm = 0;
  if(z_ij) {
    mm = k*p - n*log(exp(p) + sum_exp_p);
  }
  mm -= 0.5*precision*(p - mean)*(p - mean);
  return mm;
}


double cond_mv_mean(const NumericVector x, const NumericVector mean, const NumericMatrix Q, const int i) {
  double mm = mean[i];
  double Qii = Q(i, i);
  for(int j = 0; j < x.size(); j++) {
    if(i != j) {
      mm -= Q(i, j)*(x[j] - mean[j]) / Qii;
    }
  }
  return mm;
}

double cond_mv_mean_multinom(const NumericVector x, const IntegerVector refidx,
                             const NumericVector mean, const NumericMatrix Q,
                             const int ij) {
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

