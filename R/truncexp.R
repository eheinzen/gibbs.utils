
#' @rdname distributions
#' @export
rtruncexp <- function(n, rate = 1, a = 0, b = Inf) {
  u <- stats::runif(n)
  qtruncexp(p = u, rate = rate, a = a, b = b)
}

#' @rdname distributions
#' @export
dtruncexp <- function(x, rate = 1, a = 0, b = Inf, log = FALSE) {
  rate <- check_one_or_all(rate, length(x))
  a <- check_one_or_all(a, length(x))
  b <- check_one_or_all(b, length(x))

  flip <- rate < 0
  if(any(flip)) {
    rate[flip] <- -rate[flip]
    tmp <- a[flip]
    a[flip] <- -b[flip]
    b[flip] <- -tmp
    x[flip] <- -x[flip]
  }
  stopifnot(b > a, "An invalid bound is infinite" = is.finite(a))

  ll <- log(rate) + rate*(a - x) - log(1 - exp(-rate*(b - a)))
  r0 <- rate == 0
  ll[r0] <- -log(b[r0] - a[r0])
  ll[x < a] <- -Inf
  ll[x > b] <- -Inf

  if(log) ll else exp(ll)
}

#' @rdname distributions
#' @export
ptruncexp <- function(q, rate = 1, a = 0, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  rate <- check_one_or_all(rate, length(q))
  a <- check_one_or_all(a, length(q))
  b <- check_one_or_all(b, length(q))

  flip <- rate < 0
  if(any(flip)) {
    rate[flip] <- -rate[flip]
    tmp <- a[flip]
    a[flip] <- -b[flip]
    b[flip] <- -tmp
    q[flip] <- -q[flip]
  }
  stopifnot(b > a, "An invalid bound is infinite" = is.finite(a)) # the is.finite(a) is the only real reason I need to flip here

  num <- exp(-rate*(q - a)) - 1
  den <- exp(-rate*(b - a)) - 1

  p <- num/den
  r0 <- rate == 0
  p[r0] <- (q[r0] - a[r0])/(b[r0] - a[r0])
  p[q > b] <- 1
  p[q < a] <- 0

  if(any(flip)) {
    p[flip] <- 1 - p[flip]
  }
  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  p
}

#' @rdname distributions
#' @export
qtruncexp <- function(p, rate = 1, a = 0, b = Inf) {
  rate <- check_one_or_all(rate, length(p))
  a <- check_one_or_all(a, length(p))
  b <- check_one_or_all(b, length(p))

  flip <- rate < 0
  if(any(flip)) {
    rate[flip] <- -rate[flip]
    tmp <- a[flip]
    a[flip] <- -b[flip]
    b[flip] <- -tmp
    p[flip] <- 1 - p[flip]
  }
  stopifnot(b > a, "An invalid bound is infinite" = is.finite(a), 0 <= p, p <= 1)

  qq <- a - log(p*exp(-rate*(b - a)) - p + 1)/rate
  if(any(flip)) {
    qq[flip] <- -qq[flip]
  }
  r0 <- rate == 0
  qq[r0] <- p[r0]*(b[r0] - a[r0]) + a[r0]
  qq
}
