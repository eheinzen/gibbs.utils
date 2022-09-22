#' Sample a Poisson regression rate
#'
#' @param L the previous iteration of the log-rate
#' @param k the realized value from the Poisson distribution
#' @param mean the prior mean
#' @param precision the prior precision
#' @param ... Other arguments (not used)
#' @param nexpand The maximum number of expansions (to the right and left each)
#' @param ncontract The maximum number of contractions. If this is exceeded, the original value is returned
#' @param method The method to use to propose new values. \code{"normal"} proposes normal deviates from the current
#'   value. \code{"uniform"} proposes uniform deviates. \code{"gamma"} does moment-matching on the mean and variance
#'   to propose a conjugate gamma proposal. \code{"quadratic taylor"} proposes using a second-order Taylor
#'   approximation of the log-density at the current value, which amounts to a normal proposal with mean (usually) not equal
#'   to the current value.
#' @param width For \code{"normal"} proposals, the standard deviation(s) of proposals. For \code{"uniform"} proposals,
#' the half-width of the uniform proposal interval. For \code{"slice"}, the width of each expansion (to the right and left each).
#' For \code{"gamma"} a scaling factor to increase the variance of the proposal.
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
sample_pois_reg <- function(L, k, mean, precision, method = c("slice", "normal", "uniform", "gamma", "quadratic taylor", "mv quadratic taylor"), ...,
                            width = 1, nexpand = 10, ncontract = 100) {
  method <- match.arg(method)

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
      mean[use_norm, , drop = FALSE] + spam::rmvnorm.prec(sum(use_norm), Q = precision)
    } else matrix()
  } else {
    precision <- check_one_or_all(precision, length(L))
  }

  if(method == "slice") {
    if(is.matrix(precision)) {
      out <- slice_sample_pois_mv(L = L, k = k, k_na = is.na(k), mean = mean,
                                  Q = precision, use_norm = use_norm, norm = norm,
                                  w = width, nexpand = nexpand, ncontract = ncontract)
    } else {
      precision <- check_one_or_all(precision, length(L))
      out <- slice_sample_pois(L = L, k = k, k_na = is.na(k), mean = mean,
                               precision = precision, w = width, nexpand = nexpand, ncontract = ncontract)
    }
  } else if(method %in% c("normal", "uniform", "gamma")) {
    width <- check_one_or_all(width, length(L))
    if(method %in% c("normal", "uniform")) {
      prop <- L + if(method == "normal") stats::rnorm(length(L), 0, width) else stats::runif(length(L), -width, width)
      m <- 0L
    } else {
      prop <- width
      m <- 2L
    }

    if(is.matrix(precision)) {
      dim(prop) <- dim(L)
      out <- mh_pois_mv(method = m, L = L, proposal = prop, k = k, k_na = is.na(k), mean = mean,
                        Q = precision, use_norm = use_norm, norm = norm)
    } else {
      out <- mh_pois(method = m, L = L, proposal = prop, k = k, k_na = is.na(k), mean = mean, precision = precision)
    }
  } else if(method %in% c("quadratic taylor", "mv quadratic taylor")) {
    if(!missing(width)) warning("'width' is being ignored for this method.")

    if(is.matrix(precision)) {
      if(method == "mv quadratic taylor") {
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
        out <- mh_pois_mv(method = 1L, L = L, proposal = matrix(NA_real_), k = k, k_na = is.na(k), mean = mean,
                          Q = precision, use_norm = use_norm, norm = norm)
      }
    } else {
      if(method == "mv quadratic taylor") {
        warning("'mv quadratic taylor' is being interpreted as 'quadratic taylor' because 'precision' is not a matrix.")
      }
      out <- mh_pois(method = 1L, L = L, proposal = NA_real_, k = k, k_na = is.na(k), mean = mean, precision = precision)
    }
  }

  dim(out) <- d # could be NULL
  if(method != "slice") dim(attr(out, "accept")) <- d # could be NULL
  out
}

