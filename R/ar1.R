#' Effecient multiplication in the AR(1) setting
#'
#' @param X,Xt A matrix and its transpose
#' @param Pflat a flattened version of an AR(1) matrix. More details on this in a future version
#' @return The result of \code{t(X) \%*\% P}.
#' @export
times_flat_ar1 <- function(X, Pflat, Xt = t(X)) {
  stopifnot(dim(Pflat) == c(3, ncol(Xt)))
  times_flat_ar1_cpp(Xt, Pflat)
}
