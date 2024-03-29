#' Sample a Poisson regression rate
#'
#' @param L the previous iteration of the log-rate
#' @param k the realized value from the Poisson distribution
#' @param mean the prior mean. If this is negative infinity, \code{precision} must not
#'   be a matrix, and \code{k} must be 0.
#' @param precision the prior precision
#' @param ... Other arguments (not used)
#' @param nexpand The maximum number of expansions (to the right and left each)
#' @param ncontract The maximum number of contractions. If this is exceeded, the original value is returned
#' @param method The method to use to propose new values.
#'   \itemize{
#'     \item{\code{"normal"} proposes normal deviates from the current value.}
#'     \item{\code{"uniform"} proposes uniform deviates.}
#'     \item{\code{"gamma"} (only for Poisson) does moment-matching on the mean and variance
#'   to propose a conjugate gamma proposal.}
#'     \item{\code{"mv gamma"} (only for Poisson) does the same, but ignores off-diagonal elements of the precision matrix
#'   in the proposal (for speed's sake; this may result in low acceptance rates when elements are highly correlated).}
#'     \item{\code{"mv beta"} (only for binomial, and still somewhat experimental) does moment matching on a mean and variance
#'      approximated by the log-normal distribution (instead of logit-normal, for which there is no closed form)
#'      to propose a conjugate beta proposal; here again we ignore off-diagonal elements of the precision matrix for the proposal.}
#'     \item{\code{"quadratic taylor"} proposes using a second-order Taylor approximation of the log-density at the current value,
#'      which amounts to a normal proposal with mean (usually) not equal to the current value.}
#'      \item{\code{"mv quadratic taylor"} does the same, but uses a multivariate normal approximation.}
#'      \item{\code{"mv ind quadratic taylor"} proposes using a similar Taylor approximation, but approximates the log-density around
#'      the mean, instead of the current value. Furthermore, like \code{"mv beta"} and \code{"mv gamma"}, it ignores off-diagonal
#'      elements of the precision matrix for the proposal, which again might yield low acceptance rates when elements are highly correlated;
#'      in particular, it is not recommended for \code{sample_multinom_reg}, whose data are a priori correlated (unless you increase the
#'      variance a lot, in which case it works okay).
#'      Both of these simplifications yield a significant speed boost.}
#'      \item{\code{"mv truncated exponential"} proposes using a first-order Taylor approximation of the log-density at the current value,
#'      which amounts to a truncated exponential proposal.}
#'   }
#'   Note that \code{"mv gamma"}, \code{"mv beta"} and \code{"mv [ind ]quadratic taylor"} accept or reject an entire row at a time.
#' @param width For \code{"normal"} proposals, the standard deviation(s) of proposals. For \code{"uniform"} proposals,
#' the half-width of the uniform proposal interval. For \code{"slice"}, the width of each expansion (to the right and left each).
#' For \code{"gamma"}, \code{"beta"}, and \code{"mv ind quadratic taylor"} a scaling factor to increase the variance of the proposal
#' (it's squared, so that \code{width} is on the sd scale).
#' @param truncate Either \code{NULL} for no truncation (default), or else a list with components \code{at} to indicate the value, and \code{allow}
#'   to indicate whether values "above" (strict) or "below" (inclusive) are allowed. For convenience, negative values of \code{at} are
#'   interpreted as no truncation.
#' @param acceptance What should be the criteria for acceptance? "MH" indicates the usual Metropolis-Hastings update. "LL only" ignores the proposal
#'   densities but considers the log-likelihoods. "regardless" accepts no matter what. This is useful for testing, or for
#'   when the method is a gamma, beta, or quadratic approximation, which can be hard to accept if the initial starting point is low-density.
#' @details
#'   This function samples \code{L} conditional on \code{k}, \code{mean}, and \code{precision},
#'   where \code{k ~ Pois(exp(L))} and \code{L ~ N(mean, precision)}.
#'
#'   In the case that \code{k} is \code{NA} (akin to \code{\link{sample_binom_reg}} in the case \code{n == k == 0}),
#'   sampling is ignored in favor of a normal draw.
#'   In the special case when an entire (multivariate) row of \code{k} is \code{NA}, the entire row
#'   is simultaneously drawn (in R); when only some elements are \code{NA}, they're drawn univariately (in C++).
#'   Note that only the latter is used when \code{proposal} indicates multivariate MH.
#'
#'   This is vectorized over \code{L}, \code{k}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{L} is assumed to be multivariately distributed, and a different function is used.
#'
#'   The internals are defined in C++.
#' @seealso \url{https://en.wikipedia.org/wiki/Slice_sampling},
#' \url{https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm}, \url{https://arxiv.org/pdf/1308.0657.pdf}
#' @export
sample_pois_reg <- function(L, k, mean, precision,
                            method = c("slice", "normal", "uniform", "gamma", "mv gamma", "quadratic taylor", "mv quadratic taylor", "mv ind quadratic taylor", "mv truncated exponential"),
                            ..., width = 1, nexpand = 10, ncontract = 100, truncate = NULL, acceptance = c("MH", "LL only", "regardless")) {
  method <- match.arg(method)
  acceptance <- match.arg(acceptance)
  acceptance <- match(acceptance, c("MH", "LL only", "regardless")) - 1L

  if(length(L) != length(k)) stop("'L' and 'k' must have the same length")
  mean <- check_one_or_all(mean, length(L))
  d <- dim(L)

  if(!is.null(truncate)) {
    if(!(method %in% c("slice", "normal", "uniform"))) stop("Truncation is only supported with method = 'slice', 'normal', or 'uniform'")
    lower.tail <- truncate$allow == "below"
    at <- truncate$at
  } else {
    lower.tail <- TRUE
    at <- -1
  }
  lower.tail <- check_one_or_all(lower.tail, length(L))
  at <- check_one_or_all(at, length(L))

  if(is.matrix(precision)) {
    tmp <- is.matrix(L) + is.matrix(k)
    if(!(tmp %in% c(0, 2))) stop("Both of 'L' and 'k' must be matrices, or neither must be.")
    if(!is.matrix(L)) {
      L <- matrix(L, nrow = 1)
      k <- matrix(k, nrow = 1)
    }
    dim(lower.tail) <- dim(at) <- dim(mean) <- dim(L)
    if(any(is.infinite(mean))) stop("Infinite means aren't implemented for multivariate priors.")

    use_norm <- rowSums(!is.na(k)) == 0 & rowSums(at >= 0) == 0
    norm <- if(any(use_norm)) {
      mean[use_norm, , drop = FALSE] + spam::rmvnorm.prec(sum(use_norm), Q = precision)
    } else matrix()

  } else {
    if(method == "mv quadratic taylor") {
      warning("'mv quadratic taylor' is being interpreted as 'quadratic taylor' because 'precision' is not a matrix.")
      method <- "quadratic taylor"
    } else if(method == "mv gamma") {
      warning("'mv gamma' is being interpreted as 'gamma' because 'precision' is not a matrix.")
      method <- "gamma"
    } else if(method %in% c("mv ind quadratic taylor", "mv truncated exponential")) {
      stop(paste0("'", method, "' requires 'precision' to be a matrix"))
    }
    precision <- check_one_or_all(precision, length(L))
    inf <- is.infinite(mean)
    if(any(inf & mean > 0)) {
      stop("'mean' is positive infinite")
    } else if(any(inf)) {
      if(!all(k[inf] %in% c(0, NA_real_))) stop("If mean = -Inf, then k should be 0.")
    }
  }

  if(method == "slice") {
    if(is.matrix(precision)) {
      out <- slice_sample_pois_mv(L = L, k = k, k_na = is.na(k), mean = mean,
                                  Q = precision, use_norm = use_norm, norm = norm,
                                  trunc_at = at, lower = lower.tail,
                                  w = width, nexpand = nexpand, ncontract = ncontract)
    } else {
      precision <- check_one_or_all(precision, length(L))
      out <- slice_sample_pois(L = L, k = k, k_na = is.na(k), mean = mean, mean_inf = inf, precision = precision,
                               trunc_at = at, lower = lower.tail,
                               w = width, nexpand = nexpand, ncontract = ncontract)
    }
  } else if(method %in% c("normal", "uniform", "gamma", "quadratic taylor")) {
    if(method %in% c("normal", "uniform")) {
      width <- check_one_or_all(width, length(L))
      prop <- L + if(method == "normal") stats::rnorm(length(L), 0, width) else stats::runif(length(L), -width, width)
      m <- 0L
    } else if(method == "gamma") {
      width <- check_one_or_all(width, length(L))
      prop <- width
      m <- 2L
    } else if(method == "quadratic taylor") {
      if(!missing(width)) warning("'width' is being ignored for this method.")
      prop <- L # this is ignored
      m <- 1L
    }

    if(is.matrix(precision)) {
      dim(prop) <- dim(L)
      out <- mh_pois_mv(method = m, L = L, proposal = prop, k = k, k_na = is.na(k), mean = mean, Q = precision,
                        trunc_at = at, lower = lower.tail, use_norm = use_norm, norm = norm, acceptance = acceptance)
    } else {
      out <- mh_pois(method = m, L = L, proposal = prop, k = k, k_na = is.na(k), mean = mean, mean_inf = inf, precision = precision,
                     trunc_at = at, lower = lower.tail, acceptance = acceptance)
    }
  } else if(method %in% c("mv gamma", "mv ind quadratic taylor", "mv truncated exponential")) {
    width <- check_one_or_all(width, length(L))
    dim(width) <- dim(L)
    out <- L
    out[use_norm, ] <- norm

    not_norm <- !use_norm
    if(any(not_norm)) {
      FUN <- if(method == "mv gamma") mvgamma_pois else if(method == "mv ind quadratic taylor") mviqt_pois else mvexp_pois
      tmp <- FUN(L[not_norm, , drop = FALSE], width[not_norm, , drop = FALSE], k[not_norm, , drop = FALSE], mean[not_norm, , drop = FALSE],
                 Q = precision, acceptance = acceptance)
      out[not_norm, ] <- tmp$L
      attr(out, "accept") <- array(replace(use_norm, not_norm, tmp$accept), dim = dim(L))
    } else {
      attr(out, "accept") <- array(use_norm, dim = dim(L))
    }

  } else if(method == "mv quadratic taylor") {
    if(!missing(width)) warning("'width' is being ignored for this method.")
    use_norm[use_norm] <- seq_len(sum(use_norm))
    tmp <- t(vapply(seq_len(nrow(L)), function(i) {
      if(use_norm[i] > 0) {
        return(c(1, norm[use_norm[i], ]))
      }
      mvqt_pois(L[i, ], k[i, ], mean[i, ], precision, acceptance = acceptance)
    }, numeric(ncol(L)+1)))
    out <- tmp[, -1]
    attr(out, "accept") <- array(as.logical(tmp[, 1]), dim = dim(L))

  }

  dim(out) <- d # could be NULL
  if(method != "slice") dim(attr(out, "accept")) <- d # could be NULL
  out
}

