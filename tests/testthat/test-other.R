test_that("dmvnorm", {
  expect_equal(
    dmvnorm(1:5, 0:4, diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:5, 0:4, sqrt(100:104), log = TRUE))
  )
  expect_equal(
    dmvnorm(1:5, NULL, diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:5, 0, sqrt(100:104), log = TRUE))
  )
})

test_that("dmvlnorm", {
  expect_equal(
    dmvlnorm(1:5, 0:4, diag(1/(100:104)), log = TRUE),
    sum(dlnorm(1:5, 0:4, sqrt(100:104), log = TRUE))
  )
  expect_equal(
    dmvlnorm(1:5, NULL, diag(1/(100:104)), log = TRUE),
    sum(dlnorm(1:5, 0, sqrt(100:104), log = TRUE))
  )
})


test_that("dmatnorm", {
  expect_equal(
    dmatnorm(matrix(1:10, nrow = 5), matrix(0:9, nrow = 5), V = diag(1, 2), U = diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:10, 0:9, sqrt(c(100:104, 100:104)), log = TRUE))
  )
  expect_equal(
    dmatnorm(matrix(1:10, nrow = 5), matrix(0:9, nrow = 5), V = diag(c(1, 1/2)), U = diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:10, 0:9, sqrt(c(100:104, 2*(100:104))), log = TRUE))
  )
  expect_equal(
    dmatnorm(matrix(1:10, nrow = 5), NULL, V = diag(c(1, 1/2)), U = diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:10, 0, sqrt(c(100:104, 2*(100:104))), log = TRUE))
  )
})

test_that("*_diff", {

  Q <- matrix(c(1, -0.5, -0.5, 1.25), 2, 2)
  set.seed(88)
  x <- matrix(rnorm(20), 10, 2)
  y <- matrix(rnorm(20), 10, 2)
  mu <- matrix(rnorm(20), 10, 2)
  expect_equal(
    dmvnorm(x, mu, Q, log = TRUE) - dmvnorm(y, mu, Q, log = TRUE),
    dmvnorm_diff(x, y, mu, Q, log = TRUE)
  )

  expect_equal(
    dmvlnorm(exp(x), mu, Q, log = TRUE) - dmvlnorm(exp(y), mu, Q, log = TRUE),
    dmvlnorm_diff(exp(x), exp(y), mu, Q, log = TRUE)
  )

  tau <- 1
  set.seed(99)
  x <- exp(rnorm(10))
  y <- exp(rnorm(10))
  mu <- rnorm(10)
  expect_equal(
    dlnorm(x, meanlog = mu, sdlog = 1/sqrt(tau), log = TRUE) - dlnorm(y, meanlog = mu, sdlog = 1/sqrt(tau), log = TRUE),
    dlnorm_diff(x, y, mu, tau, log = TRUE, byrow = TRUE)
  )

  expect_equal(
    dnorm(log(x), mu, sd = 1/sqrt(tau), log = TRUE) - dnorm(log(y), mu, sd = 1/sqrt(tau), log = TRUE),
    dnorm_diff(log(x), log(y), mu, tau, log = TRUE, byrow = TRUE)
  )


  V <- matrix(c(1, -0.5, -0.5, 1.25), 2, 2)
  U <- diag(1/(100:104))
  x <- matrix(1:10, nrow = 5)
  y <- matrix(2:11, nrow = 5)
  mu <- matrix(0:9, nrow = 5)
  expect_equal(
    dmatnorm(x, mu = mu, V = V, U = U, log = TRUE) - dmatnorm(y, mu = mu, V = V, U = U, log = TRUE),
    dmatnorm_diff(x, y, mu, V = V, U = U, log = TRUE)
  )

})


test_that("*truncexp basic functionality works", {
  expect_equal(
    dtruncexp(1:10, rate = 2.5),
    stats::dexp(1:10, rate = 2.5)
  )
  expect_equal(
    dtruncexp(1:10, rate = 2.5, log = TRUE),
    stats::dexp(1:10, rate = 2.5, log = TRUE)
  )
  expect_equal(
    dtruncexp(1:2, rate = 2.5, a = 1, b = 2, log = TRUE),
    dtruncexp(2:1, rate = -2.5, a = 1, b = 2, log = TRUE)
  )
  expect_equal(
    dtruncexp(1:2, rate = 2.5, a = 1, b = 2, log = TRUE),
    dtruncexp(-(1:2), rate = -2.5, a = -2, b = -1, log = TRUE)
  )
  expect_equal(
    dtruncexp(c(0, 3), rate = 1, a = 1, b = 2),
    c(0, 0)
  )

  expect_equal(
    qtruncexp(c(0, 1), rate = 2.5),
    c(0, Inf)
  )
  expect_equal(
    qtruncexp(c(0, 1, 0, 1), rate = c(2.5, 2.5, -2.5, -2.5), a = 1, b = 2),
    c(1, 2, 1, 2)
  )

  expect_equal(
    ptruncexp(c(0, Inf), rate = 2.5),
    c(0, 1)
  )
  expect_equal(
    ptruncexp(rep(0:3, times = 2), rate = rep(c(2.5, -2.5), each = 4), a = 1, b = 2),
    rep(0:1, each = 2, times = 2)
  )
})

test_that("*truncexp errors work", {
  expect_error(dtruncexp(1, rate = -1), "An invalid bound is infinite")
  expect_error(rtruncexp(1, rate = -1), "An invalid bound is infinite")
  expect_error(qtruncexp(1, rate = -1), "An invalid bound is infinite")
  expect_error(ptruncexp(1, rate = -1), "An invalid bound is infinite")
})

test_that("rate == 0 works for *truncexp", {

  expect_equal(
    dtruncexp(c(-1, 0, 0.5, 1, 2), rate = 0, a = 0, b = 1),
    c(0, 1, 1, 1, 0)
  )
  expect_equal(
    ptruncexp(c(0, 1, 1.5, 2, 3), rate = 0, a = 1, b = 2),
    c(0, 0, 0.5, 1, 1)
  )
  expect_equal(
    qtruncexp(c(0, 0.5, 1), rate = 0, a = 1, b = 2),
    c(1, 1.5, 2)
  )
})
