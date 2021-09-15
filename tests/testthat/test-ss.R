prec <- matrix(c(10, 1, 2,
                 1, 20, 3,
                 2, 3, 30), byrow = TRUE, ncol = 3)

test_that("univariate slice sampling works for Poisson", {
  L <- 1:3
  k <- c(3, 9, 27)

  set.seed(20210729)
  one <- ss_pois_reg(L, k, mean = 2, precision = 1)
  set.seed(20210729)
  two <- ss_pois_reg(L, k, mean = rep_len(2, 3), precision = 1)
  set.seed(20210729)
  three <- ss_pois_reg(L, k, mean = 2, precision = rep_len(1, 3))
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for Poisson", {
  L <- 1:3
  k <- c(3, 9, 27)
  set.seed(20210729)
  one <- ss_pois_reg(L, k, mean = 2, precision = prec)
  set.seed(20210729)
  two <- ss_pois_reg(L, k, mean = rep_len(2, 3), precision = prec)
  set.seed(20210729)
  three <- ss_pois_reg(L, k, mean = 2, precision = prec)
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for Poisson (matrix)", {
  k <- c(3, 9, 27)
  L <- matrix(rep(1:3, each = 3), nrow = 3)
  expect_error(ss_pois_reg(L, k, mean = 2, precision = prec), "'L' and 'k' must have the same length")
  k <- matrix(rep(c(3, 9, 27), each = 3), nrow = 3)
  set.seed(20210729)
  one <- ss_pois_reg(L, k, mean = 2, precision = prec)
  set.seed(20210729)
  expect_error(ss_pois_reg(L, k, mean = rep_len(2, 3), precision = prec), "mean must be of length")
  two <- ss_pois_reg(L, k, mean = matrix(2, nrow = 3, ncol = 3), precision = prec)
  expect_equal(one, two)
  expect_true(is.matrix(one))
})






# binom -------------------------------------------------------------------


test_that("unviariate slice sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- 10
  expect_error(ss_binom_reg(p, k, n, mean = 2, precision = 1), "'p' and 'k' and 'n' must all have")
  n <- rep(10, 3)
  set.seed(20210729)
  one <- ss_binom_reg(p, k, n, mean = 2, precision = 1)
  set.seed(20210729)
  two <- ss_binom_reg(p, k, n, mean = rep_len(2, 3), precision = 1)
  set.seed(20210729)
  three <- ss_binom_reg(p, k, n, mean = 2, precision = rep_len(1, 3))
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- rep(10, 3)
  set.seed(20210729)
  one <- ss_binom_reg(p, k, n, mean = 2, precision = prec)
  set.seed(20210729)
  two <- ss_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec)
  set.seed(20210729)
  three <- ss_binom_reg(p, k, n, mean = 2, precision = prec)
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for binomial (matrix)", {
  k <- c(3, 9, 27)
  p <- matrix(rep(1:3, each = 3), nrow = 3)
  n <- matrix(10, nrow = 3, ncol = 3)
  expect_error(ss_binom_reg(p, k, n, mean = 2, precision = prec), "'p' and 'k' and 'n' must all have")
  k <- matrix(rep(c(3, 9, 27), each = 3), nrow = 3)
  set.seed(20210729)
  one <- ss_binom_reg(p, k, n, mean = 2, precision = prec)
  set.seed(20210729)
  expect_error(ss_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec), "mean must be of length")
  two <- ss_binom_reg(p, k, n, mean = matrix(2, nrow = 3, ncol = 3), precision = prec)
  expect_equal(one, two)
  expect_true(is.matrix(one))
})



# multinom -------------------------------------------------------------------

test_that("multivariate slice sampling works for binomial (matrix)", {
  tmp <- cbind(0, matrix(c(0, 0, 0.7, 0.7), nrow = 2))
  p <- array(rep(tmp, each = 400), dim = c(400, 2, 3))
  mean <- round(p)
  set.seed(99)
  n <- rpois(400*2, 1000)
  k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
  dim(k) <- c(400, 2, 3)
  stopifnot(rowSums(k, dims = 2) == n)

  z <- matrix(1, nrow = 2, ncol = 3)
  Q <- array(0, c(2, 2, 3))
  for(j in 1:3) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
  m <- apply(ss_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q), 3, mean)
  expect_true(all(abs(m - c(0, 0, log(2))) < 0.01))

})
