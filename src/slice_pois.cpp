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
double one_pois_slice_mv(NumericVector L, double k, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract) {
  double y0 = pois_LL_mv(L, k, mean, Q, i) - R::rexp(1.0);
  double left = L[i] - w;
  double right = L[i] + w;
  int j = 0;
  while(pois_LL_mv(replace_it(L, i, left), k, mean, Q, i) > y0 && j++ < nexpand) {
    left -= w;
  }
  j = 0;
  while(pois_LL_mv(replace_it(L, i, right), k, mean, Q, i) > y0 && j++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  j = 0;
  while(pois_LL_mv(replace_it(L, i, newx), k, mean, Q, i) < y0 && j++ < ncontract) {
    if(newx < L[i]) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(j == ncontract) return L[i];
  return newx;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_pois_mv(NumericMatrix L, NumericMatrix k, NumericMatrix mean, NumericMatrix Q, double w, int nexpand, int ncontract) {
  NumericMatrix out = clone(L);
  for(int r=0; r < L.nrow(); r++) {
    NumericVector kk = k(r, _);
    NumericVector mm = mean(r, _);
    for(int i=0; i < L.ncol(); i++) {
      out(r, i) = one_pois_slice_mv(out(r, _), kk[i], mm, Q, i, w, nexpand, ncontract);
    }
  }
  return out;
}
