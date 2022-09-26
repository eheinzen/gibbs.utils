
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


mvqt_pois <- function(L, k, mean, Q, accept_regardless) {

  pois_LL_mv <- function(L, k) {
    sum((k * L - exp(L))[!is.na(k)])
  }

  prop_params <- mvqt_pois_approx(L, k, mean, Q)
  proposal <- spam::rmvnorm.prec(1, mu = prop_params$mu, Q = prop_params$Q)[1, ]

  orig_params <- mvqt_pois_approx(proposal, k, mean, Q)

  if(accept_regardless) {
    accept <- TRUE
  } else {
    ratio <- pois_LL_mv(proposal, k) - pois_LL_mv(L, k) + dmvnorm_diff(proposal, L, mean, Q, log = TRUE)
    ratio <- ratio - dmvnorm(proposal, prop_params$mu, Q = prop_params$Q, log = TRUE)
    ratio <- ratio + dmvnorm(L, orig_params$mu, Q = orig_params$Q, log = TRUE)
    accept <- accept_reject(ratio)
  }
  if(accept) c(1, proposal) else c(0, L)
}


mvgamma_pois_approx <- function(k, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)
  invtau <- matrix(1.0/diag(Q), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)
  m <- exp(mean + 0.5*invtau)
  mm <- m*m
  v <- (exp(invtau) - 1.0) * mm
  alpha <- mm / v + replace(k, is.na(k), 0)
  beta <- m / v + (!is.na(k))
  gu_params(alpha = alpha, beta = beta)
}

mvgamma_pois <- function(L, mult, k, mean, Q, accept_regardless) {

  pois_LL_mv <- function(L, k) {
    rowSums(k * L - exp(L), na.rm = TRUE)
  }

  tmp <- mvgamma_pois_approx(k, mean, Q)
  alpha2 <- tmp$alpha / mult
  scale2 <- mult / tmp$beta

  proposal <- matrix(stats::rgamma(length(L), alpha2, scale = scale2), nrow = nrow(L), ncol = ncol(L))
  lproposal <- log(proposal)

  if(accept_regardless) {
    accept <- rep_len(TRUE, nrow(L))
  } else {
    ratio <- pois_LL_mv(lproposal, k) - pois_LL_mv(L, k) +
      -0.5*rowSums((lproposal + L - 2*mean) * ((lproposal - L) %*% Q))
    ratio <- ratio - rowSums((alpha2 - 1.0)*lproposal - proposal/scale2)
    ratio <- ratio + rowSums((alpha2 - 1.0)*L - exp(L)/scale2)
    accept <- vapply(ratio, accept_reject, NA)
  }

  L[accept, ] <- lproposal[accept, ]
  list(
    L = L,
    accept = accept
  )
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


mvqt_binom <- function(p, k, n, mean, Q, accept_regardless) {

  binom_LL_mv <- function(p, k, n) {
    sum(k*p - n*log(1 + exp(p)))
  }

  prop_params <- mvqt_binom_approx(p, k, n, mean, Q)
  proposal <- spam::rmvnorm.prec(1, mu = prop_params$mu, Q = prop_params$Q)[1, ]

  orig_params <- mvqt_binom_approx(proposal, k, n, mean, Q)

  if(accept_regardless) {
    accept <- TRUE
  } else {
    ratio <- binom_LL_mv(proposal, k, n) - binom_LL_mv(p, k, n) + dmvnorm_diff(proposal, p, mean, Q, log = TRUE)
    ratio <- ratio - dmvnorm(proposal, prop_params$mu, Q = prop_params$Q, log = TRUE)
    ratio <- ratio + dmvnorm(p, orig_params$mu, Q = orig_params$Q, log = TRUE)
    accept <- accept_reject(ratio)
  }
  if(accept) c(1, proposal) else c(0, p)
}


mvbeta_binom_approx <- function(k, n, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(n)) n <- matrix(n, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)
  invtau <- matrix(1.0/diag(Q), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)

  # use the log-normal approximation for when p is small
  m <- exp(mean + 0.5*invtau)
  # if(any(m > 1)) ## this is already taken into account below
  mm <- m*m
  v <- (exp(invtau) - 1.0) * mm
  if(any(v >= m*(1-m))) {
    # I think this is equivalent to 3/2 < -mu * tau
    stop("Cannot make the beta approximation")
  }

  tmp <- m*(1-m)/v - 1

  alpha <- m*tmp + k
  beta <- (1-m)*tmp + n
  gu_params(alpha = alpha, beta = beta)
}

mvbeta_binom <- function(p, mult, k, n, mean, Q, accept_regardless) {

  binom_LL_mv <- function(p, k, n) {
    rowSums(k*p - n*log(1 + exp(p)))
  }

  tmp <- mvbeta_binom_approx(k, n, mean, Q)
  alpha2 <- tmp$alpha / mult # this doesn't *exactly* multiply the variance by mult, but close enough when alpha+beta is large
  beta2 <- tmp$beta / mult

  proposal <- matrix(stats::rbeta(length(p), alpha2, beta2), nrow = nrow(p), ncol = ncol(p))
  lproposal <- logit(proposal)

  if(accept_regardless) {
    accept <- rep_len(TRUE, nrow(p))
  } else {
    ep <- expit(p)
    ratio <- binom_LL_mv(lproposal, k, n) - binom_LL_mv(p, k, n) +
      -0.5*rowSums((lproposal + p - 2*mean) * ((lproposal - p) %*% Q))
    ratio <- ratio - rowSums((alpha2 - 1)*log(proposal) + (beta2 - 1)*log(1 - proposal))
    ratio <- ratio + rowSums((alpha2 - 1)*log(ep)       + (beta2 - 1)*log(1 - ep))
    accept <- vapply(ratio, accept_reject, NA)
  }

  p[accept, ] <- lproposal[accept, ]
  list(
    p = p,
    accept = accept
  )
}

