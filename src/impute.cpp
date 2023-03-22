#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix impute_conj_mvnorm_mu_cpp(const NumericMatrix y, const NumericMatrix mu,
                                        const LogicalMatrix impute, const NumericMatrix Q,
                                        const NumericVector mu0, const NumericVector tau0) {
  NumericMatrix out = clone(mu);
  int s = 0;
  for(int i = 0; i < out.ncol(); i++) {
    double Qii = Q(i, i);
    NumericVector Qij = Q(i, _);
    for(int r = 0; r < out.nrow(); r++) {
      if(impute(r, i)) {
        double newtau = Qii + tau0[s];
        double m = tau0[s]*mu0[s] + Qii*y(r, i);
        for(int j = 0; j < out.ncol(); j++) {
          if(i != j) {
            m -= Qij[j]*(mu(r, j) - y(r, j));
          }
        }
        out(r, i) = R::rnorm(m / newtau, 1.0 / sqrt(newtau));
        s++;
      }
    }
  }
  return out;
}
