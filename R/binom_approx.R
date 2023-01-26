
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


mvqt_binom <- function(p, k, n, mean, Q, acceptance) {

  binom_LL_mv <- function(p, k, n) {
    sum(k*p - n*log(1 + exp(p)))
  }

  prop_params <- mvqt_binom_approx(p, k, n, mean, Q)
  proposal <- spam::rmvnorm.prec(1, mu = prop_params$mu, Q = prop_params$Q)[1, ]

  if(acceptance == 2) {
    accept <- TRUE
  } else {
    ratio <- binom_LL_mv(proposal, k, n) - binom_LL_mv(p, k, n) + dmvnorm_diff(proposal, p, mean, Q, log = TRUE)
    if(acceptance == 0) {
      orig_params <- mvqt_binom_approx(proposal, k, n, mean, Q)
      ratio <- ratio - dmvnorm(proposal, prop_params$mu, Q = prop_params$Q, log = TRUE)
      ratio <- ratio + dmvnorm(p, orig_params$mu, Q = orig_params$Q, log = TRUE)
    }
    accept <- accept_reject(ratio)
  }
  if(accept) c(1, proposal) else c(0, p)
}



mviqt_binom_approx <- function(around, k, n, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(n)) n <- matrix(n, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)

  tau <- matrix(diag(Q), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)

  ep <- exp(around)
  ep1 <- 1 + ep
  ep2 <- ep1*ep1

  # -H(x)
  newtau <- tau + n*ep / ep2
  newtau.inv <- 1/newtau
  # x - H^-1(x) g(x)
  newmean <- around + newtau.inv * (k - tau*(around - mean) - n*ep/ep1)

  gu_params(mu = newmean, tau = newtau)
}


mviqt_binom <- function(p, mult, k, n, mean, Q, acceptance) {

  binom_LL_mv <- function(p, k, n) {
    rowSums(k*p - n*log(1 + exp(p)))
  }

  tmp <- mviqt_binom_approx(around = mean, k, n, mean, Q) # !!! the only way we don't need to recompute params is because we approximate around the mean
  tau2 <- tmp$tau / mult
  proposal <- matrix(stats::rnorm(length(p), tmp$mu, 1/sqrt(tau2)), nrow = nrow(p), ncol = ncol(p))

  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(p))
  } else {
    ratio <- binom_LL_mv(proposal, k, n) - binom_LL_mv(p, k, n) + dmvnorm_diff(proposal, p, mu = mean, Q = Q, log = TRUE, byrow = TRUE)
    if(acceptance == 0) {
      ratio <- ratio - -0.5*rowSums((proposal + p - 2*tmp$mu) * (proposal - p) * tau2)
    }
    accept <- vapply(ratio, accept_reject, NA)
  }

  p[accept, ] <- proposal[accept, ]
  list(
    p = p,
    accept = accept
  )
}

mvbeta_binom_approx <- function(k, n, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(n)) n <- matrix(n, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)
  invtau <- matrix(1.0/diag(Q), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)

  # use the log-normal approximation for when p is small
  m <- exp(mean + 0.5*invtau)
  # if(any(m > 1)) ## this is already taken into account below
  # mm <- m*m
  # v <- (exp(invtau) - 1.0) * mm
  # if(any(v >= m*(1-m))) {
  #   # I think this is equivalent to 3/2 < -mu * tau
  #   stop("Cannot make the beta approximation")
  # }
  #
  # tmp <- m*(1-m)/v - 1
  vv <- exp(invtau) - 1.0
  if(any(vv*m >= 1-m)) {
    # I think this is equivalent to 3/2 < -mu * tau
    stop("Cannot make the beta approximation")
  }
  tmp <- (1-m)/(vv*m) - 1

  alpha <- m*tmp + k
  beta <- (1-m)*tmp + n
  gu_params(alpha = alpha, beta = beta)
}

mvbeta_binom <- function(p, mult, k, n, mean, Q, acceptance) {

  binom_LL_mv <- function(p, k, n) {
    rowSums(k*p - n*log(1 + exp(p)))
  }

  tmp <- mvbeta_binom_approx(k, n, mean, Q)
  alpha2 <- tmp$alpha / mult # this doesn't *exactly* multiply the variance by mult, but close enough when alpha+beta is large
  beta2 <- tmp$beta / mult

  proposal <- matrix(stats::rbeta(length(p), alpha2, beta2), nrow = nrow(p), ncol = ncol(p))
  lproposal <- logit(proposal)

  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(p))
  } else {
    ep <- expit(p)
    ratio <- binom_LL_mv(lproposal, k, n) - binom_LL_mv(p, k, n) + dmvnorm_diff(lproposal, p, mu = mean, Q = Q, log = TRUE, byrow = FALSE)
    ratio <- ratio - rowSums(log(proposal*(1 - proposal)) + log(ep*(1 - ep)))
    if(acceptance == 0) {
      ratio <- ratio - rowSums((alpha2 - 1)*log(proposal) + (beta2 - 1)*log(1 - proposal))
      ratio <- ratio + rowSums((alpha2 - 1)*log(ep)       + (beta2 - 1)*log(1 - ep))
    }
    accept <- vapply(ratio, accept_reject, NA)
  }

  p[accept, ] <- lproposal[accept, ]
  list(
    p = p,
    accept = accept
  )
}
