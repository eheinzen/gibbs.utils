#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_binom_ratio(double p, double proposal, double k, double n, double mean, double precision) {
  double ratio = binom_LL(proposal, k, n, mean, precision);
  ratio -= binom_LL(p, k, n, mean, precision);
  return ratio;
}


// [[Rcpp::export]]
NumericVector qt_binom_approx(double around, double k, double n, double mean, double tau) {
  // don't calculate these more than once
  double ep = exp(around);
  double ep1 = 1.0 + ep;
  double ep2 = ep1*ep1;

  // -H(x)
  double newtau = tau + n*ep / ep2;

  // x - H^-1(x) g(x)
  double newmean = around + (k - n*ep/ep1 - tau*(around - mean))/newtau;

  return NumericVector::create(newmean, 1.0/sqrt(newtau));
}


NumericVector one_qt_binom_proposal_ratio(double p, double k, double n, double mean, double precision) {

  NumericVector prop_params = qt_binom_approx(p, k, n, mean, precision);
  double proposal = R::rnorm(prop_params[0], prop_params[1]);

  NumericVector orig_params = qt_binom_approx(proposal, k, n, mean, precision);

  double ratio = one_m_binom_ratio(p, proposal, k, n, mean, precision);
  ratio -= R::dnorm(proposal, prop_params[0], prop_params[1], 1);
  ratio += R::dnorm(p, orig_params[0], orig_params[1], 1);

  return NumericVector::create(proposal, ratio);
}


// [[Rcpp::export]]
NumericVector mh_binom(bool qt, NumericVector p, NumericVector proposal, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision) {
  NumericVector out = clone(p);
  LogicalVector accept(p.size());

  for(int i=0; i < p.size(); i++) {
    if(n[i] == 0.0) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      accept[i] = true;
      continue;
    }

    double ratio, prop;
    if(qt) { // 'proposal' is ignored
      NumericVector prop_ratio = one_qt_binom_proposal_ratio(p[i], k[i], n[i], mean[i], precision[i]);
      prop = prop_ratio[0];
      ratio = prop_ratio[1];
    } else {
      prop = proposal[i];
      ratio = one_m_binom_ratio(p[i], prop, k[i], n[i], mean[i], precision[i]);
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
NumericVector mh_binom_mv(bool qt, NumericMatrix p, NumericMatrix proposal, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q,
                         LogicalVector use_norm, NumericMatrix norm) {
  NumericMatrix out = clone(p);
  LogicalMatrix accept(p.nrow(), p.ncol());

  int nrm = 0;
  for(int r=0; r < p.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      for(int i=0; i < p.ncol(); i++) {
        accept(r, i) = true;
      }
      continue;
    }

    NumericVector kk = k(r, _);
    NumericVector nn = n(r, _);
    NumericVector mm = mean(r, _);
    for(int i = 0; i < p.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);

      if(nn[i] == 0.0) {
        out(r, i) = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        accept(r, i) = true;
        continue;
      }

      double ratio, prop;
      if(qt) { // 'proposal' is ignored
        NumericVector prop_ratio = one_qt_binom_proposal_ratio(out(r, i), kk[i], nn[i], mmm, Q(i, i));
        prop = prop_ratio[0];
        ratio = prop_ratio[1];
      } else {
        prop = proposal(r, i);
        ratio = one_m_binom_ratio(out(r, i), prop, kk[i], nn[i], mmm, Q(i, i));
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
