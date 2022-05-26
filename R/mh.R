#' @rdname pois_reg
#' @param proposal The method to use to propose new values. \code{"normal"} proposes normal deviates from the current
#'   value. \code{"uniform"} proposes uniform deviates. \code{"quadratic taylor"} proposes using a second-order Taylor
#'   approximation of the log-density at the current value, which amounts to a normal proposal with mean (usually) not equal
#'   to the current value.
#' @param proposal_sd For \code{"normal"} proposals, the standard deviation(s) of proposals. For \code{"uniform"} proposals,
#' the half-width of the uniform proposal interval.
#' @export
mh_pois_reg <- function(L, k, mean, precision, proposal = c("normal", "uniform", "quadratic taylor", "mv quadratic taylor"), ..., proposal_sd = 1) {
  proposal <- match.arg(proposal)

  if(length(L) != length(k)) stop("'L' and 'k' must have the same length")
  mean <- check_one_or_all(mean, length(L))
  d <- dim(L)

  if(is.matrix(precision)) {
    tmp <- is.matrix(L) + is.matrix(k)
    if(!(tmp %in% c(0, 2))) stop("Both of 'L' and 'k' must be matrices, or neither must be.")
    if(!is.matrix(L)) {
      L <- matrix(L, nrow = 1)
      k <- matrix(k, nrow = 1)
    }
    dim(mean) <- dim(L)

    use_norm <- rowSums(!is.na(k)) == 0
    norm <- if(any(use_norm)) {
      mean[use_norm, , drop = FALSE] + chol_mvrnorm(sum(use_norm), mu = 0, Precision = precision)
    } else matrix()
  } else {
    precision <- check_one_or_all(precision, length(L))
  }

  if(proposal %in% c("normal", "uniform")) {
    proposal_sd <- check_one_or_all(proposal_sd, length(L))
    prop <- L + if(proposal == "normal") {
      stats::rnorm(length(L), 0, proposal_sd)
    } else stats::runif(length(L), -proposal_sd, proposal_sd)

    if(is.matrix(precision)) {
      dim(prop) <- dim(L)
      out <- m_pois_mv(L, prop, k, is.na(k), mean, precision, use_norm = use_norm, norm = norm)
    } else {
      out <- m_pois(L, prop, k, is.na(k), mean, precision)
    }
  } else if(proposal %in% c("quadratic taylor", "mv quadratic taylor")) {
    if(!missing(proposal_sd)) warning("'proposal_sd' is being ignored for this method.")

    if(is.matrix(precision)) {
      if(proposal == "mv quadratic taylor") {
        use_norm[use_norm] <- seq_len(sum(use_norm))
        tmp <- t(vapply(seq_len(nrow(L)), function(i) {
          if(use_norm[i] > 0) {
            return(c(1, norm[use_norm[i], ]))
          }
          mvqt_pois(L[i, ], k[i, ], mean[i, ], precision)
        }, numeric(ncol(L)+1)))
        out <- tmp[, -1]
        attr(out, "accept") <- array(as.logical(tmp[, 1]), dim = dim(L))

      } else {
        out <- qt_pois_mv(L, k, is.na(k), mean, precision, use_norm = use_norm, norm = norm)
      }
    } else {
      if(proposal == "mv quadratic taylor") {
        warning("'mv quadratic taylor' is being interpreted as 'quadratic taylor' because 'precision' is not a matrix.")
      }
      out <- qt_pois(L, k, is.na(k), mean, precision)
    }
  }

  dim(attr(out, "accept")) <- dim(out) <- d # could be NULL
  out
}




mvqt_pois_approx <- function(around, k, mean, Q) {
  ep <- exp(around)

  # -H(x)
  newQ <- Q + diag(replace(ep, is.na(k), 0), length(around))
  newQ.inv <- chol_inv(newQ)
  # x - H^-1(x) g(x)
  newmean <- around + newQ.inv %*% (replace(k - ep, is.na(k), 0) - Q %*% (around - mean))

  list(newmean, newQ, newQ.inv)
}


mvqt_pois <- function(L, k, mean, Q) {

  pois_LL_mv <- function(L, k, mean, Q) {
    sum((k * L - exp(L))[!is.na(k)]) - 0.5*as.vector(t(L - mean) %*% Q %*% (L - mean))
  }

  prop_params <- mvqt_pois_approx(L, k, mean, Q);
  proposal <- chol_mvrnorm(1, prop_params[[1]], Sigma = prop_params[[3]]);

  orig_params <- mvqt_pois_approx(proposal, k, mean, Q);

  ratio <- pois_LL_mv(proposal, k, mean, Q)
  ratio <- ratio - dmvnorm(proposal, prop_params[[1]], Q = prop_params[[2]], log = TRUE)

  ratio <- ratio - pois_LL_mv(L, k, mean, Q)
  ratio <- ratio + dmvnorm(L, orig_params[[1]], Q = orig_params[[2]], log = TRUE)
  accept <- accept_reject(ratio)
  # print(exp(ratio))
  # list(if(accept) proposal else L, accept, exp(ratio), proposal,  prop_params, orig_params);
  if(accept) c(1, proposal) else c(0, L)
}







#' @rdname binom_reg
#' @export
mh_binom_reg <- function(p, k, n, mean, precision, proposal = c("normal", "uniform", "quadratic taylor", "mv quadratic taylor"), ..., proposal_sd = 1) {
  proposal <- match.arg(proposal)

  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  stopifnot(k <= n)
  mean <- check_one_or_all(mean, length(p))
  d <- dim(p)

  if(is.matrix(precision)) {
    tmp <- is.matrix(p) + is.matrix(k) + is.matrix(n)
    if(!(tmp %in% c(0, 3))) stop("All of 'p', 'k', and 'n' must be matrices, or none must be.")
    if(!is.matrix(p)) {
      p <- matrix(p, nrow = 1)
      k <- matrix(k, nrow = 1)
      n <- matrix(n, nrow = 1)
    }
    dim(mean) <- dim(p)

    use_norm <- rowSums(n) == 0
    norm <- if(any(use_norm)) {
      mean[use_norm, , drop = FALSE] + chol_mvrnorm(sum(use_norm), mu = 0, Precision = precision)
    } else matrix()
  } else {
    precision <- check_one_or_all(precision, length(p))
  }

  if(proposal %in% c("normal", "uniform")) {
    proposal_sd <- check_one_or_all(proposal_sd, length(p))
    prop <- p + if(proposal == "normal") {
      stats::rnorm(length(p), 0, proposal_sd)
    } else stats::runif(length(p), -proposal_sd, proposal_sd)

    if(is.matrix(precision)) {
      dim(prop) <- dim(p)
      out <- m_binom_mv(p, prop, k, n, mean, precision, use_norm = use_norm, norm = norm)
    } else {
      out <- m_binom(p, prop, k, n, mean, precision)
    }
  } else if(proposal %in% c("quadratic taylor", "mv quadratic taylor")) {
    if(!missing(proposal_sd)) warning("'proposal_sd' is being ignored for this method.")

    if(is.matrix(precision)) {
      if(proposal == "mv quadratic taylor") {
        use_norm[use_norm] <- seq_len(sum(use_norm))
        tmp <- t(vapply(seq_len(nrow(p)), function(i) {
          if(use_norm[i] > 0) {
            return(c(1, norm[use_norm[i], ]))
          }
          mvqt_binom(p[i, ], k[i, ], n[i, ], mean[i, ], precision)
        }, numeric(ncol(p)+1)))
        out <- tmp[, -1]
        attr(out, "accept") <- array(as.logical(tmp[, 1]), dim = dim(p))

      } else {
        out <- qt_binom_mv(p, k, n, mean, precision, use_norm = use_norm, norm = norm)
      }
    } else {
      if(proposal == "mv quadratic taylor") {
        warning("'mv quadratic taylor' is being interpreted as 'quadratic taylor' because 'precision' is not a matrix.")
      }
      out <- qt_binom(p, k, n, mean, precision)
    }
  }

  dim(attr(out, "accept")) <- dim(out) <- d # could be NULL
  out
}



mvqt_binom_approx <- function(around, k, n, mean, Q) {
  ep <- exp(around)
  ep1 <- 1 + ep
  ep2 <- ep1*ep1

  # -H(x)
  newQ <- Q + diag(n*ep / ep2, length(around))
  newQ.inv <- chol_inv(newQ)
  # x - H^-1(x) g(x)
  newmean <- around + newQ.inv %*% (k - Q %*% (around - mean) - n*ep/ep1)

  list(newmean, newQ, newQ.inv)
}


mvqt_binom <- function(p, k, n, mean, Q) {

  binom_LL_mv <- function(p, k, n, mean, Q) {
    sum(k*p - n*log(1 + exp(p))) - 0.5*as.vector(t(p - mean) %*% Q %*% (p - mean))
  }

  prop_params <- mvqt_binom_approx(p, k, n, mean, Q);
  proposal <- chol_mvrnorm(1, prop_params[[1]], Sigma = prop_params[[3]]);

  orig_params <- mvqt_binom_approx(proposal, k, n, mean, Q);

  ratio <- binom_LL_mv(proposal, k, n, mean, Q)
  ratio <- ratio - dmvnorm(proposal, prop_params[[1]], Q = prop_params[[2]], log = TRUE)

  ratio <- ratio - binom_LL_mv(p, k, n, mean, Q)
  ratio <- ratio + dmvnorm(p, orig_params[[1]], Q = orig_params[[2]], log = TRUE)
  accept <- accept_reject(ratio)
  # print(exp(ratio))
  # list(if(accept) proposal else p, accept, exp(ratio), proposal,  prop_params, orig_params);
  if(accept) c(1, proposal) else c(0, p)
}

