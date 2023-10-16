#' Impute Conjugate Values
#'
#' @param y A matrix of values drawn from the multivariate normal distribution
#'   whose mean we are to impute.
#' @param mu the means we wish to impute. Should be the same size as \code{y}.
#'   No NA's are allowed; indices to impute are passed to \code{impute}. Instead,
#'   usually values are initialized and allowed to update within the MCMC framework.
#' @param impute A logical matrix the size of \code{y}, indicating where to impute.
#' @param Q the precision of the multivariate normal distribution from which \code{y} comes.
#'   Should be square, with dimension \code{ncol(y)}.
#' @param mu0,tau0 The prior means/precisions for \code{mu}.
#'  Can be a matrix the same size as \code{mu}, in which case
#'  only the values corresponding to \code{impute == TRUE} are used, or a vector with length
#'  \code{sum(impute)}.
#' @name impute
NULL

#' @rdname impute
#' @export
impute_conj_mvnorm_mu <- function(y, mu, impute, Q, mu0, tau0) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  if(!is.matrix(mu)) mu <- matrix(mu, nrow = 1)
  if(!is.matrix(impute)) impute <- matrix(impute, nrow = 1)
  if(!identical(dim(y), dim(mu)) || !identical(dim(mu), dim(impute)))
    stop("dim(mu) != dim(y) or dim(mu) != dim(impute)")

  N <- sum(impute)
  mu0 <- if(is.matrix(mu0)) mu0[impute] else check_one_or_all(mu0, N)
  tau0 <- if(is.matrix(tau0)) tau0[impute] else check_one_or_all(tau0, N)

  impute_conj_mvnorm_mu_cpp(y = y, mu = mu, impute = impute, Q = Q, mu0 = mu0, tau0 = tau0)
}
