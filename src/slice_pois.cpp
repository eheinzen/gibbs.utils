#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_pois_slice(double L, double k, double mean, double precision, double w, int nexpand, int ncontract) {
  double y0 = pois_LL(L, k, mean, precision) - R::rexp(1.0);
  double left = L - w;
  double right = L + w;
  int i = 0;
  while(pois_LL(left, k, mean, precision) > y0 && i++ < nexpand) {
    left -= w;
  }
  i = 0;
  while(pois_LL(right, k, mean, precision) > y0 && i++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  i = 0;
  while(pois_LL(newx, k, mean, precision) < y0 && i++ < ncontract) {
    if(newx < L) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(i == ncontract) return L;
  return newx;
}


// [[Rcpp::export]]
NumericVector slice_sample_pois(NumericVector L, NumericVector k, NumericVector mean, NumericVector precision, double w, int nexpand, int ncontract) {
  NumericVector out(L.size());
  for(int i=0; i < L.size(); i++) {
    out[i] = one_pois_slice(L[i], k[i], mean[i], precision[i], w, nexpand, ncontract);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_pois_mv(NumericMatrix L, NumericMatrix k, NumericMatrix mean, NumericMatrix Q, double w, int nexpand, int ncontract) {
  NumericMatrix out = clone(L);
  for(int r=0; r < L.nrow(); r++) {
    NumericVector kk = k(r, _);
    NumericVector mm = mean(r, _);
    for(int i=0; i < L.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);
      out(r, i) = one_pois_slice(out(r, i), kk[i], mmm, Q(i, i), w, nexpand, ncontract);
    }
  }
  return out;
}
