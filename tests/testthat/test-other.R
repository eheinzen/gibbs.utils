test_that("dmvnorm", {
  expect_equal(
    dmvnorm(1:5, 0:4, diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:5, 0:4, sqrt(100:104), log = TRUE))
  )

  Q <- matrix(c(1, -0.5, -0.5, 1.25), 2, 2)
  set.seed(88)
  x <- matrix(rnorm(20), 10, 2)
  y <- matrix(rnorm(20), 10, 2)
  mu <- matrix(rnorm(20), 10, 2)
  expect_equal(
    dmvnorm(x, mu, Q, log = TRUE) - dmvnorm(y, mu, Q, log = TRUE),
    dmvnorm_diff(x, y, mu, Q, log = TRUE)
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
})
