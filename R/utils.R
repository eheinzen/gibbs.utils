
#' @useDynLib gibbs.utils
#' @importFrom Rcpp sourceCpp
NULL

check_one_or_all <- function(x, len) {
  nm <- deparse(substitute(x))
  if(length(x) == 1) rep_len(x, len) else if(length(x) == len) x else stop(nm, " must be of length 1 or ", len)
}

#' Other Utilities
#'
#' @param x A matrix for \code{chol_inv} or a numeric vector for \code{expit}
#' @param p A numeric vector
#' @rdname utilities
#' @export
chol_inv <- function(x) {
  chol2inv(chol(x))
}

#' @rdname utilities
#' @export
logit <- function(p) {
  log(p) - log(1 - p)
}

#' @rdname utilities
#' @export
expit <- function(x) {
  1 / (1 + exp(-x))
}

