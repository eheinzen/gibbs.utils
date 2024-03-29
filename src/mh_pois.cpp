#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_pois_ratio(const double L, const double proposal, const double k,
                        const double mean, const double precision,
                        const double trunc_at, const bool lower,
                        const int acceptance) {
  if(acceptance == 2) return 0.0;
  return pois_LL(proposal, k, mean, precision, trunc_at, lower) - pois_LL(L, k, mean, precision, trunc_at, lower);
}


void qt_pois_approx(const double around, const double k,
                    const double mean, const double tau,
                    double& outmean, double& outsd) {
  // don't calculate this more than once
  double ep = exp(around);

  // -H(x)
  double newtau = tau + ep;
  outsd = 1.0/sqrt(newtau);

  // x - H^-1(x) g(x)
  double newmean = around + (k - ep - tau*(around - mean))/newtau;
  outmean = newmean;
}


void one_qt_pois_proposal_ratio(const double L, const double k,
                                const double mean, const double precision,
                                double& outproposal, double& outratio,
                                const int acceptance) {

  double prop_mean, prop_sd;
  qt_pois_approx(L, k, mean, precision, prop_mean, prop_sd);
  double proposal = R::rnorm(prop_mean, prop_sd);
  outproposal = proposal;

  double ratio = one_m_pois_ratio(L, proposal, k, mean, precision, -1, true, acceptance);
  if(acceptance == 0) {
    double orig_mean, orig_sd;
    qt_pois_approx(proposal, k, mean, precision, orig_mean, orig_sd);
    ratio -= R::dnorm(proposal, prop_mean, prop_sd, 1);
    ratio += R::dnorm(L, orig_mean, orig_sd, 1);
  }
  outratio = ratio;
}

void gamma_pois_approx(const double k, const double mean, const double tau,
                       double& outalpha, double& outbeta) {
  double invtau = 1.0/tau;
  double m = exp(mean + 0.5*invtau);
  //double mm = m*m;
  //double v = (exp(invtau) - 1.0) * mm;
  //double alpha = mm / v + k;
  //double beta = m / v + 1.0;
  double vv = exp(invtau) - 1.0;
  double alpha = 1.0/vv + k;
  double beta = 1.0/(m*vv) + 1.0;
  outalpha = alpha;
  outbeta = beta;
}

void one_gamma_pois_proposal_ratio(const double L, const double mult, const double k,
                                   const double mean, const double precision,
                                   double& outproposal, double& outratio,
                                   const int acceptance) {

  double alpha, beta;
  gamma_pois_approx(k, mean, precision, alpha, beta);
  double alpha2 = alpha/mult, scale2 = mult/beta;
  double proposal = R::rgamma(alpha2, scale2);
  double lproposal = log(proposal);
  outproposal = lproposal;

  double ratio = one_m_pois_ratio(L, lproposal, k, mean, precision, -1, true, acceptance) - (lproposal - L);
  if(acceptance == 0) {
    ratio -= (alpha2 - 1.0)*lproposal - proposal/scale2;
    ratio += (alpha2 - 1.0)*L - exp(L)/scale2;
  }
  outratio = ratio;
}



// [[Rcpp::export]]
NumericVector mh_pois(const int method, const NumericVector L, const NumericVector proposal,
                      const NumericVector k, const LogicalVector k_na,
                      const NumericVector mean, const LogicalVector mean_inf, const NumericVector precision,
                      const NumericVector trunc_at, const LogicalVector lower,
                      const int acceptance) {
  NumericVector out = clone(L);
  LogicalVector accept(L.size());

  for(int i=0; i < L.size(); i++) {
    if(k_na[i]) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      accept[i] = true;
      continue;
    } else if(mean_inf[i]) {
      out[i] = R_NegInf;
      accept[i] = true;
      continue;
    }

    double ratio, prop;
    if(method == 1) { // 'proposal' is ignored
      if(k_na[i]) stop("Something went wrong");
      one_qt_pois_proposal_ratio(L[i], k[i], mean[i], precision[i], prop, ratio, acceptance);
    } else if(method == 2) {
      if(k_na[i]) stop("Something went wrong");
      one_gamma_pois_proposal_ratio(L[i], proposal[i], k[i], mean[i], precision[i], prop, ratio, acceptance);
    } else {
      prop = proposal[i];
      ratio = one_m_pois_ratio(L[i], prop, k[i], mean[i], precision[i], trunc_at[i], lower[i], acceptance);
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
NumericVector mh_pois_mv(const int method, const NumericMatrix L, const NumericMatrix proposal,
                         const NumericMatrix k, const LogicalMatrix k_na,
                         const NumericMatrix mean, const NumericMatrix Q,
                         const NumericMatrix trunc_at, const LogicalMatrix lower,
                         const LogicalVector use_norm, const NumericMatrix norm,
                         const int acceptance) {
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
    NumericVector aa = trunc_at(r, _);
    LogicalVector bb = lower(r, _);
    NumericVector tmp = out(r, _);
    NumericVector propp = proposal(r, _);
    LogicalVector acc = accept(r, _);
    for(int i = 0; i < L.ncol(); i++) {
      // safe not to replace tmp[i] because it's not used in this calculation
      double mmm = cond_mv_mean(tmp, mm, Q, i);

      if(kk_na[i]) {
        tmp[i] = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        acc[i] = true;
        continue;
      }

      double ratio, prop;
      if(method == 1) { // 'proposal' is ignored
        if(kk_na[i]) stop("Something went wrong");
        one_qt_pois_proposal_ratio(tmp[i], kk[i], mmm, Q(i, i), prop, ratio, acceptance);
      } else if(method == 2) {
        if(kk_na[i]) stop("Something went wrong");
        one_gamma_pois_proposal_ratio(tmp[i], propp[i], kk[i], mmm, Q(i, i), prop, ratio, acceptance);
      } else {
        prop = propp[i];
        ratio = one_m_pois_ratio(tmp[i], prop, kk[i], mmm, Q(i, i), aa[i], bb[i], acceptance);
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
