#include <Rcpp.h>
using namespace Rcpp;


double binom_LL(double p, double k, double n, double mean, double precision) {
  double pp = log(p);
  double qq = log(1.0 - p);
  return k*pp + (n-k)*qq - pp - qq - 0.5*precision*(pp - qq - mean)*(pp - qq - mean);
}

double expit(double x) {
  return 1.0 / (1.0 + exp(-x));
}
double logit(double p) {
  return log(p) - log(1.0 - p);
}

// [[Rcpp::export]]
double one_binom_slice(double p, double k, double n, double mean, double precision, double w, int nexpand, int ncontract) {
  double y0 = binom_LL(p, k, n, mean, precision) - R::rexp(1.0);
  double left = logit(p) - w;
  double right = logit(p) + w;
  int i = 0;
  while(binom_LL(expit(left), k, n, mean, precision) > y0 && i++ < ncontract) {
    left -= w;
  }
  i = 0;
  while(binom_LL(expit(right), k, n, mean, precision) > y0 && i++ < ncontract) {
    right += w;
  }
  left = expit(left);
  right = expit(right);
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

