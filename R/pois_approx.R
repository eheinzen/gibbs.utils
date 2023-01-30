
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


mvqt_pois <- function(L, k, mean, Q, acceptance) {

  pois_LL_mv <- function(L, k) {
    sum((k * L - exp(L))[!is.na(k)])
  }

  prop_params <- mvqt_pois_approx(L, k, mean, Q)
  proposal <- spam::rmvnorm.prec(1, mu = prop_params$mu, Q = prop_params$Q)[1, ]


  if(acceptance == 2) {
    accept <- TRUE
  } else {
    ratio <- pois_LL_mv(proposal, k) - pois_LL_mv(L, k) + dmvnorm_diff(proposal, L, mean, Q, log = TRUE)
    if(acceptance == 0) {
      orig_params <- mvqt_pois_approx(proposal, k, mean, Q)
      ratio <- ratio - dmvnorm(proposal, prop_params$mu, Q = prop_params$Q, log = TRUE)
      ratio <- ratio + dmvnorm(L, orig_params$mu, Q = orig_params$Q, log = TRUE)
    }
    accept <- accept_reject(ratio)
  }
  if(accept) c(1, proposal) else c(0, L)
}


mviqt_pois_approx <- function(around, k, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)

  tau <- matrix(diag(Q), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)

  ep <- exp(around)

  # -H(x)
  newtau <- tau + replace(ep, is.na(k), 0)
  # x - H^-1(x) g(x)
  newmean <- around + (replace(k - ep, is.na(k), 0) - tau*(around - mean))/newtau

  gu_params(mu = newmean, tau = newtau)
}


mviqt_pois <- function(L, mult, k, mean, Q, acceptance) {

  pois_LL_mv <- function(L, k) {
    rowSums(k * L - exp(L), na.rm = TRUE)
  }

  tmp <- mviqt_pois_approx(around = mean, k, mean, Q) # !!! the only way we don't need to recompute params is because we approximate around the mean
  tau2 <- tmp$tau / mult^2
  proposal <- matrix(stats::rnorm(length(L), tmp$mu, 1/sqrt(tau2)), nrow = nrow(L), ncol = ncol(L))
  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(L))
  } else {
    ratio <- pois_LL_mv(proposal, k) - pois_LL_mv(L, k) + dmvnorm_diff(proposal, L, mu = mean, Q = Q, log = TRUE, byrow = TRUE)
    if(acceptance == 0) {
      ratio <- ratio - -0.5*rowSums((proposal + L - 2*tmp$mu) * (proposal - L) * tau2)
    }
    accept <- vapply(ratio, accept_reject, NA)
  }

  L[accept, ] <- proposal[accept, ]
  list(
    L = L,
    accept = accept
  )
}


mvgamma_pois_approx <- function(k, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)
  invtau <- matrix(1.0/diag(Q), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)
  m <- exp(mean + 0.5*invtau)
  # mm <- m*m
  # v <- (exp(invtau) - 1.0) * mm
  # alpha <- mm / v + replace(k, is.na(k), 0)
  # beta <- m / v + (!is.na(k))
  vv <- exp(invtau) - 1.0
  alpha <- 1/vv + replace(k, is.na(k), 0)
  beta <- 1/(m*vv) + (!is.na(k))
  if(any(!is.finite(alpha)) || any(!is.finite(beta))) stop("Can't seem to make the mv gamma approximation. Is your precision too small?")
  gu_params(alpha = alpha, beta = beta)
}

mvgamma_pois <- function(L, mult, k, mean, Q, acceptance) {

  pois_LL_mv <- function(L, k) {
    rowSums(k * L - exp(L), na.rm = TRUE)
  }

  tmp <- mvgamma_pois_approx(k, mean, Q)
  alpha2 <- tmp$alpha / mult^2
  scale2 <- mult^2 / tmp$beta

  proposal <- matrix(stats::rgamma(length(L), shape = alpha2, scale = scale2), nrow = nrow(L), ncol = ncol(L))
  lproposal <- log(proposal)

  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(L))
  } else {
    eL <- exp(L)
    ratio <- pois_LL_mv(lproposal, k) - pois_LL_mv(L, k) + dmvlnorm_diff(proposal, eL, mu = mean, Q = Q, log = TRUE, byrow = TRUE)
    if(acceptance == 0) {
      ratio <- ratio - rowSums((alpha2 - 1.0)*lproposal - proposal/scale2)
      ratio <- ratio + rowSums((alpha2 - 1.0)*L - exp(L)/scale2)
    }
    accept <- vapply(ratio, accept_reject, NA)
  }

  L[accept, ] <- lproposal[accept, ]
  list(
    L = L,
    accept = accept
  )
}




mvexp_pois_approx <- function(around, k, mean, Q) {
  if(!is.matrix(k)) k <- matrix(k, nrow = 1)
  if(!is.matrix(mean)) mean <- matrix(mean, nrow = 1)

  ep <- exp(around)

  # g(x)
  g <- replace(k - ep, is.na(k), 0) - (around - mean) %*% Q

  gu_params(slope = g)
}


mvexp_pois <- function(L, mult, k, mean, Q, acceptance) {

  pois_LL_mv <- function(L, k) {
    rowSums(k * L - exp(L), na.rm = TRUE)
  }

  tmp <- mvexp_pois_approx(L, k, mean, Q)

  aL <- L - mult
  bL <- L + mult
  proposal <- matrix(
    rtruncexp(
      n = nrow(L)*ncol(L),
      rate = tmp$slope,
      a = aL,
      b = bL
    ),
    nrow = nrow(L), ncol = ncol(L)
  )
  if(acceptance == 2) {
    accept <- rep_len(TRUE, nrow(L))
  } else {
    ratio <- pois_LL_mv(proposal, k) - pois_LL_mv(L, k) + dmvnorm_diff(proposal, L, mu = mean, Q = Q, log = TRUE, byrow = TRUE)
    if(acceptance == 0) {
      tmp2 <- mvexp_pois_approx(proposal, k, mean, Q)
      aprop <- proposal - mult
      bprop <- proposal + mult
      ratio <- ratio -
        rowSums(dtruncexp(proposal, rate = tmp$slope, a = aL, b = bL, log = TRUE)) +
        rowSums(dtruncexp(L, rate = tmp2$slope, a = aprop, b = bprop, log = TRUE))
    }
    accept <- vapply(ratio, accept_reject, NA)
  }

  L[accept, ] <- proposal[accept, ]
  list(
    L = L,
    accept = accept
  )
}

