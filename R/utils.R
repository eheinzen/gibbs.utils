
#' @useDynLib gibbs.utils
#' @importFrom Rcpp sourceCpp
NULL

check_one_or_all <- function(x, len) {
  nm <- deparse(substitute(x))
  if(length(x) == 1) rep_len(x, len) else if(length(x) == len) x else stop(nm, " must be of length 1 or ", len)
}

my_inv <- function(x) {
  chol2inv(chol(x))
}


logit <- function(p) log(p) - log(1 - p)
expit <- function(x) 1 / (1 + exp(-x))

