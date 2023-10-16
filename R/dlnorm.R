
npint <- function(x) x < 0 | x != floor(x)

#' @rdname distributions
#' @export
rdlnorm <- function(n, meanlog = 0, sdlog = 1) {
  floor(stats::rlnorm(n = n, meanlog = meanlog, sdlog = sdlog))
}

#' @rdname distributions
#' @export
ddlnorm <- function(x, meanlog = 0, sdlog = 1, log = FALSE) {
  idx <- npint(x)
  if(any(idx)) warning("'x' contains non-positive-integers")
  out <- stats::plnorm(x + 1, meanlog = meanlog, sdlog = sdlog) -
    stats::plnorm(x, meanlog = meanlog, sdlog = sdlog)
  out[idx] <- 0
  if(log) log(out) else out
}

#' @rdname distributions
#' @export
pdlnorm <- function(x, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
  stats::plnorm(floor(x) + 1, meanlog = meanlog, sdlog = sdlog, lower.tail = lower.tail, log.p = log.p)
}

