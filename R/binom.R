#' Sample a binomial regression rate
#'
#' @param p the previous iteration of the logit-probability
#' @param k the realized value from the binomial distribution
#' @param n the number of trials
#' @inheritParams sample_pois_reg
#' @inherit sample_pois_reg seealso
#' @details
#'   This function samples \code{p} conditional on \code{k}, \code{n}, \code{mean}, and \code{precision},
#'   where \code{k ~ Binom(n, expit(p))} and \code{p ~ N(mean, precision)}.
#'
#'   In the case that \code{n} is zero, sampling is ignored in favor of a normal draw.
#'   In the special case when an entire (multivariate) row of \code{n} is zero, the entire row
#'   is simultaneously drawn (in R); when only some elements are zero, they're drawn univariately (in C++).
#'   Note that only the latter is used when \code{proposal} indicates multivariate MH.
#'
#'   This is vectorized over \code{p}, \code{k}, \code{n}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{p} is assumed to be multivariately distributed, and different internals are used.
#'
#'   The internals are defined in C++.
#' @export
sample_binom_reg <- function(p, k, n, mean, precision,
                             method = c("slice", "normal", "uniform", "quadratic taylor", "mv quadratic taylor", "mv ind quadratic taylor", "mv beta"),
                             ..., width = 1, nexpand = 10, ncontract = 100, acceptance = c("MH", "LL only", "regardless")) {
  method <- match.arg(method)
  acceptance <- match.arg(acceptance)
  acceptance <- match(acceptance, c("MH", "LL only", "regardless")) - 1L

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
      mean[use_norm, , drop = FALSE] + spam::rmvnorm.prec(sum(use_norm), Q = precision)
    } else matrix()
  } else {
    if(method == "mv quadratic taylor") {
      warning("'mv quadratic taylor' is being interpreted as 'quadratic taylor' because 'precision' is not a matrix.")
      method <- "quadratic taylor"
    } else if(method %in% c("mv beta", "mv ind quadratic taylor")) {
      stop(paste0("'", method, "' requires 'precision' to be a matrix"))
    }
    precision <- check_one_or_all(precision, length(p))
  }


  if(method == "slice") {
    if(is.matrix(precision)) {
      out <- slice_sample_binom_mv(p = p, k = k, n = n, mean = mean,
                                   Q = precision, use_norm = use_norm, norm = norm,
                                   w = width, nexpand = nexpand, ncontract = ncontract)
    } else {
      out <- slice_sample_binom(p = p, k = k, n = n, mean = mean,
                                precision = precision, w = width, nexpand = nexpand, ncontract = ncontract)
    }
  } else if(method %in% c("normal", "uniform", "quadratic taylor")) {
    if(method %in% c("normal", "uniform")) {
      width <- check_one_or_all(width, length(p))
      prop <- p + if(method == "normal") stats::rnorm(length(p), 0, width) else stats::runif(length(p), -width, width)
      qt <- FALSE
    } else if(method == "quadratic taylor") {
      if(!missing(width)) warning("'width' is being ignored for this method.")
      prop <- p # this is ignored
      qt <- TRUE
    }
    if(is.matrix(precision)) {
      dim(prop) <- dim(p)
      out <- mh_binom_mv(qt = qt, p = p, proposal = prop, k = k, n = n, mean = mean,
                         Q = precision, use_norm = use_norm, norm = norm, acceptance = acceptance)
    } else {
      out <- mh_binom(qt = qt, p = p, proposal = prop, k = k, n = n, mean = mean, precision = precision, acceptance = acceptance)
    }
  } else if(method %in% c("mv beta", "mv ind quadratic taylor")) {
    width <- check_one_or_all(width, length(p))
    dim(width) <- dim(p)
    out <- p
    out[use_norm, ] <- norm

    not_norm <- !use_norm
    if(any(not_norm)) {
      FUN <- if(method == "mv beta") mvbeta_binom else mviqt_binom
      tmp <- FUN(p[not_norm, , drop = FALSE], width[not_norm, , drop = FALSE], k[not_norm, , drop = FALSE], n[not_norm, , drop = FALSE],
                 mean[not_norm, , drop = FALSE], Q = precision, acceptance = acceptance)
      out[not_norm, ] <- tmp$p
      attr(out, "accept") <- array(replace(use_norm, not_norm, tmp$accept), dim = dim(p))
    } else {
      attr(out, "accept") <- array(use_norm, dim = dim(p))
    }

  } else if(method == "mv quadratic taylor") {
    if(!missing(width)) warning("'width' is being ignored for this method.")
    use_norm[use_norm] <- seq_len(sum(use_norm))
    tmp <- t(vapply(seq_len(nrow(p)), function(i) {
      if(use_norm[i] > 0) {
        return(c(1, norm[use_norm[i], ]))
      }
      mvqt_binom(p[i, ], k[i, ], n[i, ], mean[i, ], precision, acceptance = acceptance)
    }, numeric(ncol(p)+1)))
    out <- tmp[, -1]
    attr(out, "accept") <- array(as.logical(tmp[, 1]), dim = dim(p))
  }

  dim(out) <- d # could be NULL
  if(method != "slice") dim(attr(out, "accept")) <- d # could be NULL
  out
}
