
#' @rdname binom_reg
#' @param proposal_sd The standard deviation(s) of proposals
#' @export
mh_binom_reg <- function(p, k, n, mean, precision, ..., proposal_sd = 1) {
  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  mean <- check_one_or_all(mean, length(p))
  proposal_sd <- check_one_or_all(proposal_sd, length(p))
  d <- dim(p)

  if(is.matrix(precision)) {
    FUN <- mh_binom_mv
    tmp <- is.matrix(p) + is.matrix(k) + is.matrix(n)
    if(!(tmp %in% c(0, 3))) stop("All of 'p', 'k', and 'n' must be matrices, or none must be.")
    if(!is.matrix(p)) {
      p <- matrix(p, nrow = 1)
      k <- matrix(k, nrow = 1)
      n <- matrix(n, nrow = 1)
    }
    proposal <- rnorm(length(p), p, proposal_sd)
    dim(proposal) <- dim(mean) <- dim(p)
  } else {
    precision <- check_one_or_all(precision, length(p))
    FUN <- mh_binom
    proposal <- rnorm(length(p), p, proposal_sd)
  }
  out <- FUN(p, proposal, k, n, mean, precision)
  dim(out) <- d # could be NULL
  out
}
