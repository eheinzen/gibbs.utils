
#' Draw Gibbs Samples
#'
#' Sample from the posterior conditional on all but one parameter (at a time).
#'
#' @inheritParams conjugacy
#' @param a0,b0 Optional limits of truncation for a truncated normal prior.
#' @param ... Passed through to \code{\link{conj_matlm_beta}}
#' @name gibbs
NULL

#' @rdname gibbs
#' @export
gs_matlm_beta <- function(beta, ...,
                              a0 = rep_len(-Inf, length(beta)), b0 = rep_len(Inf, length(beta)),
                              params.only = FALSE) {
  params <- conj_matlm_beta(..., params.only = TRUE)
  if(params.only) {
    tmp <- lapply(seq_along(beta), function(i) {
      cond_mvnorm(y = beta, mu = params$mu, Q = params$Q, which = i, a = a0[i], b = b0[i], params.only = TRUE)
    })
    return(gu_params(mu = lapply(tmp, "[[", "mu"), sd = lapply(tmp, "[[", "sd"), a = a0, b = b0))
  }
  out <- beta
  for(i in seq_along(out)) {
    out[i] <- cond_mvnorm(y = out, mu = params$mu, Q = params$Q, which = i, a = a0[i], b = b0[i])
  }
  out
}
