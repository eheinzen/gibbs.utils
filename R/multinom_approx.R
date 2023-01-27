
mviqt_multinom_approx <- function(around, z, which_i, is_ref, k, n, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(n)) n <- matrix(n, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)

  tau <- matrix(diag(Q), nrow = nrow(k), ncol = ncol(Q), byrow = TRUE)

  ep <- exp(around) * z

  sum_ep <- vapply(seq_len(max(which_i)), function(i) rowSums(ep[, which_i == i, drop = FALSE]), numeric(nrow(k)))
  dim(sum_ep) <- c(nrow(k), max(which_i))
  sum_ep <- sum_ep[, which_i, drop = FALSE]
  ep1 <- sum_ep
  sum_ep <- sum_ep - ep
  ep2 <- ep1*ep1

  # -H(x)
  newtau <- tau + (n*ep*sum_ep / ep2)[, !is_ref]

  # x - H^-1(x) g(x)
  newmean <- around[, !is_ref] + ((k - n*ep/ep1)[, !is_ref] - tau*(around[, !is_ref] - mean))/newtau
  gu_params(mu = newmean, tau = newtau)
}


mviqt_multinom <- function(p, mult, z, which_i, is_ref, k, n, mean, Q, acceptance) {

  multinom_LL_mv <- function(x) {
    kz <- k*z
    rowSums(kz*x) - rowSums(vapply(seq_len(max(which_i)), function(i) {
      nn <- rowSums(kz[, which_i == i, drop = FALSE])
      stopifnot(n[, which_i == i & is_ref] == nn)
      e <- log(rowSums(exp(x[, which_i == i, drop = FALSE]) * z[, which_i == i, drop = FALSE]))
      nn * e
    }, numeric(nrow(k))))
  }

  around <- p
  around[, !is_ref] <- mean
  tmp <- mviqt_multinom_approx(around, z, which_i, is_ref, k, n, mean, Q) # !!! the only way we don't need to recompute params is because we approximate around the mean
  tau2 <- tmp$tau / mult^2
  proposal <- p
  proposal[, !is_ref] <- matrix(stats::rnorm(nrow(p)*sum(!is_ref), tmp$mu, 1/sqrt(tau2)), nrow = nrow(p), ncol = sum(!is_ref))

  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(p))
  } else {
    ratio <- multinom_LL_mv(proposal) - multinom_LL_mv(p) + dmvnorm_diff(proposal[, !is_ref], p[, !is_ref], mu = mean, Q = Q, log = TRUE, byrow = TRUE)
    if(acceptance == 0) {
      ratio <- ratio - -0.5*rowSums((proposal[, !is_ref] + p[, !is_ref] - 2*tmp$mu) * (proposal[, !is_ref] - p[, !is_ref]) * tau2)
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
