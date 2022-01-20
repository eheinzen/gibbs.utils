
#' @rdname binom_reg
#' @param proposal The method to use to propose new values. \code{"normal"} proposes normal deviates from the current
#'   value. \code{"uniform"} proposes uniform deviates. \code{"quadratic taylor"} proposes using a second-order Taylor
#'   approximation of the log-density at the current value, which amounts to a normal proposal with mean (usually) not equal
#'   to the current value.
#' @seealso \url{https://arxiv.org/pdf/1308.0657.pdf}
#' @param proposal_sd For \code{"normal"} proposals, the standard deviation(s) of proposals. For \code{"uniform"} proposals,
#' the half-width of the uniform proposal interval.
#' @export
mh_binom_reg <- function(p, k, n, mean, precision, proposal = c("normal", "uniform", "quadratic taylor"), ..., proposal_sd = 1) {
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

    use_norm <- n == 0
    norm <- if(any(use_norm)) {
      mean[use_norm] + stats::rnorm(sum(use_norm), mean = 0, sd = 1 / sqrt(precision[use_norm]))
    } else numeric()
  }

  if(proposal %in% c("normal", "uniform")) {
    proposal_sd <- check_one_or_all(proposal_sd, length(p))
    prop <- p + if(proposal == "normal") {
      stats::rnorm(length(p), 0, proposal_sd)
    } else stats::runif(length(p), -proposal_sd, proposal_sd)

    if(is.matrix(precision)) {
      FUN <- m_binom_mv
      dim(prop) <- dim(p)
    } else FUN <- m_binom
    out <- FUN(p, prop, k, n, mean, precision, use_norm = use_norm, norm = norm)

  } else if(proposal == "quadratic taylor") {
    if(!missing(proposal_sd)) warning("'proposal_sd' is being ignored for this method.")

    if(is.matrix(precision)) {
      FUN <- qt_binom_mv
    } else FUN <- qt_binom
    out <- FUN(p, k, n, mean, precision, use_norm = use_norm, norm = norm)
  }

  dim(attr(out, "accept")) <- dim(out) <- d # could be NULL
  out
}
