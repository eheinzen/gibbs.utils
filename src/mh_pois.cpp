#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_pois_ratio(double L, double proposal, double k, double mean, double precision) {
  double ratio = pois_LL(proposal, k, mean, precision);
  ratio -= pois_LL(L, k, mean, precision);
  return ratio;
}


// [[Rcpp::export]]
NumericVector qt_pois_approx(double around, double k, double mean, double tau) {
  // don't calculate this more than once
  double ep = exp(around);

  // -H(x)
  double newtau = tau + ep;

  // x - H^-1(x) g(x)
  double newmean = around + (k - ep - tau*(around - mean))/newtau;

  return NumericVector::create(newmean, 1.0/sqrt(newtau));
}


NumericVector one_qt_pois_proposal_ratio(double L, double k, double mean, double precision) {

  NumericVector prop_params = qt_pois_approx(L, k, mean, precision);
  double proposal = R::rnorm(prop_params[0], prop_params[1]);

  NumericVector orig_params = qt_pois_approx(proposal, k, mean, precision);
  double ratio = one_m_pois_ratio(L, proposal, k, mean, precision);
  ratio -= R::dnorm(proposal, prop_params[0], prop_params[1], 1);
  ratio += R::dnorm(L, orig_params[0], orig_params[1], 1);

  return NumericVector::create(proposal, ratio);
}


// [[Rcpp::export]]
NumericVector mh_pois(bool qt, NumericVector L, NumericVector proposal, NumericVector k, LogicalVector k_na, NumericVector mean, NumericVector precision) {
  NumericVector out = clone(L);
  LogicalVector accept(L.size());

  for(int i=0; i < L.size(); i++) {
    if(k_na[i]) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      accept[i] = true;
      continue;
    }

    double ratio, prop;
    if(qt) { // 'proposal' is ignored
      NumericVector prop_ratio = one_qt_pois_proposal_ratio(L[i], k[i], mean[i], precision[i]);
      ratio = prop_ratio[1];
      prop = prop_ratio[0];
    } else {
      prop = proposal[i];
      ratio = one_m_pois_ratio(L[i], prop, k[i], mean[i], precision[i]);
    }
    bool a = accept_reject(ratio);
    accept[i] = a;
    if(a) {
      out[i] = prop;
    }
  }
  out.attr("accept") = accept;
  return out;
}


// [[Rcpp::export]]
NumericVector mh_pois_mv(bool qt, NumericMatrix L, NumericMatrix proposal, NumericMatrix k, LogicalMatrix k_na, NumericMatrix mean, NumericMatrix Q,
                          LogicalVector use_norm, NumericMatrix norm) {
  NumericMatrix out = clone(L);
  LogicalMatrix accept(L.nrow(), L.ncol());

  int nrm = 0;
  for(int r=0; r < L.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      for(int i=0; i < L.ncol(); i++) {
        accept(r, i) = true;
      }
      continue;
    }

    NumericVector kk = k(r, _);
    LogicalVector kk_na = k_na(r, _);
    NumericVector mm = mean(r, _);
    for(int i = 0; i < L.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);

      if(kk_na[i]) {
        out(r, i) = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        accept(r, i) = true;
        continue;
      }

      double ratio, prop;
      if(qt) { // 'proposal' is ignored
        NumericVector prop_ratio = one_qt_pois_proposal_ratio(out(r, i), kk[i], mmm, Q(i, i));
        ratio = prop_ratio[1];
        prop = prop_ratio[0];
      } else {
        prop = proposal(r, i);
        ratio = one_m_pois_ratio(out(r, i), proposal(r, i), kk[i], mmm, Q(i, i));
      }
      bool a = accept_reject(ratio);
      accept(r, i) = a;
      if(a) {
        out(r, i) = prop;
      }
    }
  }
  out.attr("accept") = accept;
  return out;
}
