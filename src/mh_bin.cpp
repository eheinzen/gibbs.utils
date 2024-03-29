#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_binom_ratio(const double p, const double proposal,
                         const double k, const double n,
                         const double mean, const double precision,
                         const int acceptance) {
  if(acceptance == 2) return 0.0;
  return binom_LL(proposal, k, n, mean, precision) - binom_LL(p, k, n, mean, precision);
}


void qt_binom_approx(const double around, const double k, const double n,
                     const double mean, const double tau,
                     double& outmean, double& outsd) {
  // don't calculate these more than once
  double ep = exp(around);
  double ep1 = 1.0 + ep;
  double ep2 = ep1*ep1;

  // -H(x)
  double newtau = tau + n*ep / ep2;
  outsd = 1.0/sqrt(newtau);

  // x - H^-1(x) g(x)
  double newmean = around + (k - n*ep/ep1 - tau*(around - mean))/newtau;
  outmean = newmean;
}

void one_qt_binom_proposal_ratio(const double p, const double k, const double n,
                                 const double mean, const double precision,
                                 double& outproposal, double& outratio,
                                 const int acceptance) {

  double prop_mean, prop_sd;
  qt_binom_approx(p, k, n, mean, precision, prop_mean, prop_sd);
  double proposal = R::rnorm(prop_mean, prop_sd);
  outproposal = proposal;

  double ratio = one_m_binom_ratio(p, proposal, k, n, mean, precision, acceptance);
  if(acceptance == 0) {
    double orig_mean, orig_sd;
    qt_binom_approx(proposal, k, n, mean, precision, orig_mean, orig_sd);
    ratio -= R::dnorm(proposal, prop_mean, prop_sd, 1);
    ratio += R::dnorm(p, orig_mean, orig_sd, 1);
  }
  outratio = ratio;
}


// [[Rcpp::export]]
NumericVector mh_binom(const bool qt, const NumericVector p, const NumericVector proposal,
                       const NumericVector k, const NumericVector n,
                       const NumericVector mean, const NumericVector precision,
                       const int acceptance) {
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
      one_qt_binom_proposal_ratio(p[i], k[i], n[i], mean[i], precision[i], prop, ratio, acceptance);
    } else {
      prop = proposal[i];
      ratio = one_m_binom_ratio(p[i], prop, k[i], n[i], mean[i], precision[i], acceptance);
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
NumericVector mh_binom_mv(const bool qt, const NumericMatrix p, const NumericMatrix proposal,
                          const NumericMatrix k, const NumericMatrix n,
                          const NumericMatrix mean, const NumericMatrix Q,
                          const LogicalVector use_norm, const NumericMatrix norm,
                          const int acceptance) {
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
    NumericVector tmp = out(r, _);
    NumericVector propp = proposal(r, _);
    LogicalVector acc = accept(r, _);
    for(int i = 0; i < p.ncol(); i++) {
      // safe not to replace tmp[i] because it's not used in this calculation
      double mmm = cond_mv_mean(tmp, mm, Q, i);

      if(nn[i] == 0.0) {
        tmp[i] = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        acc[i] = true;
        continue;
      }

      double ratio, prop;
      if(qt) { // 'proposal' is ignored
        one_qt_binom_proposal_ratio(tmp[i], kk[i], nn[i], mmm, Q(i, i), prop, ratio, acceptance);
      } else {
        prop = propp[i];
        ratio = one_m_binom_ratio(tmp[i], prop, kk[i], nn[i], mmm, Q(i, i), acceptance);
      }
      bool a = accept_reject(ratio);
      acc[i] = a;
      if(a) {
        tmp[i] = prop;
      }
    }
    out(r, _) = tmp;
    accept(r, _) = acc;
  }
  out.attr("accept") = accept;
  return out;
}
