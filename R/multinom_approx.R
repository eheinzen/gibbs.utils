
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




mvexp_multinom_approx <- function(around, z, which_i, is_ref, k, n, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(n)) n <- matrix(n, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)

  ep <- exp(around) * z

  sum_ep <- vapply(seq_len(max(which_i)), function(i) rowSums(ep[, which_i == i, drop = FALSE]), numeric(nrow(k)))
  dim(sum_ep) <- c(nrow(k), max(which_i))
  ep1 <- sum_ep[, which_i, drop = FALSE]

  # g(x)
  g <- (k - n*ep/ep1)[, !is_ref] - (around[, !is_ref] - mean) %*% Q
  gu_params(slope = g)
}

mvexp_multinom <- function(p, mult, z, which_i, is_ref, k, n, mean, Q, acceptance) {

  multinom_LL_mv <- function(x) {
    kz <- k*z
    rowSums(kz*x) - rowSums(vapply(seq_len(max(which_i)), function(i) {
      nn <- rowSums(kz[, which_i == i, drop = FALSE])
      stopifnot(n[, which_i == i & is_ref] == nn)
      e <- log(rowSums(exp(x[, which_i == i, drop = FALSE]) * z[, which_i == i, drop = FALSE]))
      nn * e
    }, numeric(nrow(k))))
  }

  tmp <- mvexp_multinom_approx(p, z, which_i, is_ref, k, n, mean, Q)

  ap <- p[, !is_ref] - mult
  bp <- p[, !is_ref] + mult

  proposal <- p
  proposal[, !is_ref] <- matrix(rtruncexp(
    n = nrow(p)*sum(!is_ref),
    rate = tmp$slope,
    a = ap,
    b = bp
  ), nrow = nrow(p), ncol = sum(!is_ref))

  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(p))
  } else {
    ratio <- multinom_LL_mv(proposal) - multinom_LL_mv(p) + dmvnorm_diff(proposal[, !is_ref], p[, !is_ref], mu = mean, Q = Q, log = TRUE, byrow = TRUE)
    if(acceptance == 0) {
      tmp2 <- mvexp_multinom_approx(proposal, z, which_i, is_ref, k, n, mean, Q)
      aprop <- proposal[, !is_ref] - mult
      bprop <- proposal[, !is_ref] + mult
      ratio <- ratio -
        rowSums(dtruncexp(proposal[, !is_ref], rate = tmp$slope, a = ap, b = bp, log = TRUE)) +
        rowSums(dtruncexp(p[, !is_ref], rate = tmp2$slope, a = aprop, b = bprop, log = TRUE))
    }
    accept <- vapply(ratio, accept_reject, NA)
  }

  p[accept, ] <- proposal[accept, ]
  list(
    p = p,
    accept = accept
  )
}















