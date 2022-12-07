#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector times_flat_ar1_cpp(NumericMatrix xt, NumericMatrix flat_ar1) {
  NumericMatrix out(xt.nrow(), xt.ncol());
  for(int i=0; i < xt.ncol(); i++) {
    out(_, i) = xt(_, i)*flat_ar1(1, i);
    if(i > 0) {
      out(_, i) = out(_, i) + xt(_, i-1)*flat_ar1(0, i);
    }
    if(i < xt.ncol()-1) {
      out(_, i) = out(_, i) + xt(_, i+1)*flat_ar1(2, i);
    }
  }
  return out;
}

