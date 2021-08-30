#' Draw Conditional Samples
#'
#' Sample from a conditional distribution.
#'
#' @param y A vector of values drawn from the multivariate normal distribution. Element \code{which} is ignored.
#' @param which Which index to get the conditional distribution of
#' @param Q the precision of the multivariate normal distribution from which \code{y} comes
#' @param mu the mean of the normal distribution from which \code{y} comes
#' @param params.only Should just a list of the updated parameters be returned?
#' @rdname conditional
#' @export
cond_mvnorm <- function(y, mu, Q, which, params.only = FALSE) {
  if(length(which) != 1) stop("This is only implemented for one parameter at a time")
  Q11 <- Q[which, which]
  mean <- mu[which] - sum(Q[-which, which]/Q11 * (y - mu)[-which])

  if(params.only) return(gu_params(mu = mean, tau = Q11, sd = 1 / sqrt(Q11)))
  stats::rnorm(1, mean = mean, sd = 1 / sqrt(Q11))
}

