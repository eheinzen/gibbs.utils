#include <Rcpp.h>
using namespace Rcpp;


double binom_LL(double p, double k, double n, double mean, double precision) {
  return -k*log(1.0 + exp(-p)) - (n-k)*log(1.0 + exp(p)) - 0.5*precision*(p - mean)*(p - mean);
}

// [[Rcpp::export]]
double one_binom_slice(double p, double k, double n, double mean, double precision, double w, int nexpand, int ncontract) {
  double y0 = binom_LL(p, k, n, mean, precision) - R::rexp(1.0);
  double left = p - w;
  double right = p + w;
  int i = 0;
  while(binom_LL(left, k, n, mean, precision) > y0 && i++ < ncontract) {
    left -= w;
  }
  i = 0;
  while(binom_LL(right, k, n, mean, precision) > y0 && i++ < ncontract) {
    right += w;
  }
  double newx = R::runif(left, right);
  i = 0;
  while(binom_LL(newx, k, n, mean, precision) < y0 && i++ < ncontract) {
    if(newx < p) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(i == ncontract) return p;
  return newx;
}




// [[Rcpp::export]]
NumericVector slice_sample_binom(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision,
                                 double w, int nexpand, int ncontract) {
  NumericVector out(p.size());
  for(int i=0; i < p.size(); i++) {
    out[i] = one_binom_slice(p[i], k[i], n[i], mean[i], precision[i], w, nexpand, ncontract);
  }
  return out;
}





double binom_LL_mv(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericMatrix Q) {
  double mm = 0.0;
  for(int i = 0; i < p.size(); i++) {
    mm += -k[i]*log(1.0 + exp(-p[i])) - (n[i]-k[i])*log(1.0 + exp(p[i])) - 0.5*sum((p - mean) * Q(_, i)) * (p[i] - mean[i]);
  }
  return mm;
}

NumericVector replace_itt(NumericVector x, int i, double value) {
  NumericVector out = clone(x);
  out[i] = value;
  return out;
}

// [[Rcpp::export]]
double one_binom_slice_mv(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericMatrix Q, int i, double w, int nexpand, int ncontract) {
  double y0 = binom_LL_mv(p, k, n, mean, Q) - R::rexp(1.0);
  double left = p[i] - w;
  double right = p[i] + w;
  int j = 0;
  while(binom_LL_mv(replace_itt(p, i, left), k, n, mean, Q) > y0 && j++ < nexpand) {
    left -= w;
  }
  j = 0;
  while(binom_LL_mv(replace_itt(p, i, right), k, n, mean, Q) > y0 && j++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  j = 0;
  while(binom_LL_mv(replace_itt(p, i, newx), k, n, mean, Q) < y0 && j++ < ncontract) {
    if(newx < p[i]) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(j == ncontract) return p[i];
  return newx;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_binom_mv(NumericMatrix p, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q, double w, int nexpand, int ncontract) {
  NumericMatrix out = clone(p);
  for(int r=0; r < p.nrow(); r++) {
    for(int i=0; i < p.ncol(); i++) {
      out(r, i) = one_binom_slice_mv(out(r, _), k(r, _), n(r, _), mean(r, _), Q, i, w, nexpand, ncontract);
    }
  }
  return out;
}
