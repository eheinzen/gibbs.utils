
#' Helper functions for approximations
#'
#' @param k,n Number of observations
#' @param around The vector around which to perform an approximation.
#' @param mean,Q Priors
#' @name approximations
NULL

#' @rdname approximations
#' @export
mvqt_pois_approx <- function(around, k, mean, Q) {
  ep <- exp(around)

  # -H(x)
  newQ <- Q + diag(replace(ep, is.na(k), 0), length(around))
  newQ.inv <- chol_inv(newQ)
  # x - H^-1(x) g(x)
  newmean <- around + newQ.inv %*% (replace(k - ep, is.na(k), 0) - Q %*% (around - mean))

  gu_params(mu = newmean, Q = newQ, Q.inv = newQ.inv)
}


mvqt_pois <- function(L, k, mean, Q) {

  pois_LL_mv <- function(L, k, mean, Q) {
    sum((k * L - exp(L))[!is.na(k)]) - 0.5*as.vector(t(L - mean) %*% Q %*% (L - mean))
  }

  prop_params <- mvqt_pois_approx(L, k, mean, Q)
  proposal <- spam::rmvnorm.prec(1, mu = prop_params$mu, Q = prop_params$Q)[1, ]

  orig_params <- mvqt_pois_approx(proposal, k, mean, Q)

  ratio <- pois_LL_mv(proposal, k, mean, Q)
  ratio <- ratio - dmvnorm(proposal, prop_params$mu, Q = prop_params$Q, log = TRUE)

  ratio <- ratio - pois_LL_mv(L, k, mean, Q)
  ratio <- ratio + dmvnorm(L, orig_params$mu, Q = orig_params$Q, log = TRUE)
  accept <- accept_reject(ratio)
  # print(exp(ratio))
  # list(if(accept) proposal else L, accept, exp(ratio), proposal,  prop_params, orig_params);
  if(accept) c(1, proposal) else c(0, L)
}


mvgamma_pois_approx <- function(k, mean, Q) {
  invtau <- 1.0/diag(Q)
  m <- exp(mean + 0.5*invtau)
  mm <- m*m
  v <- (exp(invtau) - 1.0) * mm
  alpha <- mm / v + replace(k, is.na(k), 0)
  beta <- m / v + (!is.na(k))
  gu_params(alpha = alpha, beta = beta)
}

mvgamma_pois <- function(L, mult, k, mean, Q) {

  pois_LL_mv <- function(L, k, mean, Q) {
    sum((k * L - exp(L))[!is.na(k)]) - 0.5*as.vector(t(L - mean) %*% Q %*% (L - mean))
  }

  tmp <- mvgamma_pois_approx(k, mean, Q)
  alpha2 <- tmp$alpha / mult
  scale2 <- mult / tmp$beta

  proposal <- rgamma(length(L), alpha2, scale = scale2)
  lproposal <- log(proposal)

  ratio <- pois_LL_mv(lproposal, k, mean, Q)
  ratio <- ratio - sum((alpha2 - 1.0)*lproposal - proposal/scale2)
  ratio <- ratio - pois_LL_mv(L, k, mean, Q)
  ratio <- ratio + sum((alpha2 - 1.0)*L - exp(L)/scale2)
  accept <- accept_reject(ratio)
  if(accept) c(1, lproposal) else c(0, L)
}


#' @rdname approximations
#' @export
mvqt_binom_approx <- function(around, k, n, mean, Q) {
  ep <- exp(around)
  ep1 <- 1 + ep
  ep2 <- ep1*ep1

  # -H(x)
  newQ <- Q + diag(n*ep / ep2, length(around))
  newQ.inv <- chol_inv(newQ)
  # x - H^-1(x) g(x)
  newmean <- around + newQ.inv %*% (k - Q %*% (around - mean) - n*ep/ep1)

  gu_params(mu = newmean, Q = newQ, Q.inv = newQ.inv)
}


mvqt_binom <- function(p, k, n, mean, Q) {

  binom_LL_mv <- function(p, k, n, mean, Q) {
    sum(k*p - n*log(1 + exp(p))) - 0.5*as.vector(t(p - mean) %*% Q %*% (p - mean))
  }

  prop_params <- mvqt_binom_approx(p, k, n, mean, Q)
  proposal <- spam::rmvnorm.prec(1, mu = prop_params$mu, Q = prop_params$Q)[1, ]

  orig_params <- mvqt_binom_approx(proposal, k, n, mean, Q)

  ratio <- binom_LL_mv(proposal, k, n, mean, Q)
  ratio <- ratio - dmvnorm(proposal, prop_params$mu, Q = prop_params$Q, log = TRUE)

  ratio <- ratio - binom_LL_mv(p, k, n, mean, Q)
  ratio <- ratio + dmvnorm(p, orig_params$mu, Q = orig_params$Q, log = TRUE)
  accept <- accept_reject(ratio)
  # print(exp(ratio))
  # list(if(accept) proposal else p, accept, exp(ratio), proposal,  prop_params, orig_params);
  if(accept) c(1, proposal) else c(0, p)
}

