#include <Rcpp.h>
using namespace Rcpp;


double pois_LL(double L, double k, double mean, double precision) {
  return k*L - exp(L) - 0.5*precision*(L - mean)*(L - mean);
}

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





double pois_LL_mv(NumericVector L, NumericVector k, NumericVector mean, NumericMatrix Q) {
  double mm = 0.0;
  for(int i = 0; i < L.size(); i++) {
    mm += k[i]*L[i] - exp(L[i]) - 0.5*sum((L - mean) * Q(_, i)) * (L[i] - mean[i]);
  }
  return mm;
}

NumericVector replace_it(NumericVector x, int i, double value) {
  NumericVector out = clone(x);
  out[i] = value;
  return out;
}

// [[Rcpp::export]]
double one_pois_slice_mv(NumericVector L, NumericVector k, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract) {
  double y0 = pois_LL_mv(L, k, mean, Q) - R::rexp(1.0);
  double left = L[i] - w;
  double right = L[i] + w;
  int j = 0;
  while(pois_LL_mv(replace_it(L, i, left), k, mean, Q) > y0 && j++ < nexpand) {
    left -= w;
  }
  j = 0;
  while(pois_LL_mv(replace_it(L, i, right), k, mean, Q) > y0 && j++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  j = 0;
  while(pois_LL_mv(replace_it(L, i, newx), k, mean, Q) < y0 && j++ < ncontract) {
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
NumericVector slice_sample_pois_mv(NumericVector L, NumericVector k, NumericVector mean, NumericMatrix Q, double w, int nexpand, int ncontract) {
  NumericVector out = clone(L);
  for(int i=0; i < L.size(); i++) {
    out[i] = one_pois_slice_mv(out, k, mean, Q, i, w, nexpand, ncontract);
  }
  return out;
}
