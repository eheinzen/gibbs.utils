prec <- matrix(c(10, 1, 2,
                 1, 20, 3,
                 2, 3, 30), byrow = TRUE, ncol = 3)

test_that("univariate slice sampling works for Poisson", {
  L <- 1:3
  k <- c(3, 9, 27)

  set.seed(20210729)
  one <- sample_pois_reg(L, k, mean = 2, precision = 1, method = "slice")
  set.seed(20210729)
  two <- sample_pois_reg(L, k, mean = rep_len(2, 3), precision = 1, method = "slice")
  set.seed(20210729)
  three <- sample_pois_reg(L, k, mean = 2, precision = rep_len(1, 3), method = "slice")
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for Poisson", {
  L <- 1:3
  k <- c(3, 9, 27)
  set.seed(20210729)
  one <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "slice")
  set.seed(20210729)
  two <- sample_pois_reg(L, k, mean = rep_len(2, 3), precision = prec, method = "slice")
  set.seed(20210729)
  three <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "slice")
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for Poisson (matrix)", {
  k <- c(3, 9, 27)
  L <- matrix(rep(1:3, each = 3), nrow = 3)
  expect_error(sample_pois_reg(L, k, mean = 2, precision = prec, method = "slice"), "'L' and 'k' must have the same length")
  k <- matrix(rep(c(3, 9, 27), each = 3), nrow = 3)
  set.seed(20210729)
  one <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "slice")
  set.seed(20210729)
  expect_error(sample_pois_reg(L, k, mean = rep_len(2, 3), precision = prec, method = "slice"), "'x' must be of length")
  two <- sample_pois_reg(L, k, mean = matrix(2, nrow = 3, ncol = 3), precision = prec, method = "slice")
  expect_equal(one, two)
  expect_true(is.matrix(one))
})


test_that("Poisson slice sampling works for NAs", {
  m <- L <- rep_len(0, 10)
  k <- c(rep_len(NA, 5), 1:5)
  tau <- 10
  set.seed(20211203)
  norm <- rnorm(5, 0, sd = 1/sqrt(tau))
  set.seed(20211203)
  ss <- sample_pois_reg(L, k, m, tau, method = "slice")
  expect_equal(ss[1:5], norm)

  m <- L <- matrix(1, 4, 3)
  k <- matrix(c(rep_len(NA, 6), 1:6), byrow = TRUE, nrow = 4)
  tau <- diag(10, 3)
  set.seed(20211203)
  norm <- spam::rmvnorm.prec(2, rep(1, 3), Q = tau)
  set.seed(20211203)
  ss <- sample_pois_reg(L, k, m, tau, method = "slice")
  expect_equal(ss[1:2, ], norm)
})





# binom -------------------------------------------------------------------


test_that("unviariate slice sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- 10
  expect_error(sample_binom_reg(p, k, n, mean = 2, precision = 1, method = "slice"), "'p' and 'k' and 'n' must all have")
  n <- rep(10, 3)
  set.seed(20210729)
  one <- sample_binom_reg(p, k, n, mean = 2, precision = 1, method = "slice")
  set.seed(20210729)
  two <- sample_binom_reg(p, k, n, mean = rep_len(2, 3), precision = 1, method = "slice")
  set.seed(20210729)
  three <- sample_binom_reg(p, k, n, mean = 2, precision = rep_len(1, 3), method = "slice")
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- rep(10, 3)
  set.seed(20210729)
  one <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "slice")
  set.seed(20210729)
  two <- sample_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec, method = "slice")
  set.seed(20210729)
  three <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "slice")
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate slice sampling works for binomial (matrix)", {
  k <- c(3, 9, 27)
  p <- matrix(rep(1:3, each = 3), nrow = 3)
  n <- matrix(100, nrow = 3, ncol = 3)
  expect_error(sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "slice"), "'p' and 'k' and 'n' must all have")
  k <- matrix(rep(c(3, 9, 27), each = 3), nrow = 3)
  set.seed(20210729)
  one <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "slice")
  set.seed(20210729)
  expect_error(sample_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec, method = "slice"), "'x' must be of length")
  two <- sample_binom_reg(p, k, n, mean = matrix(2, nrow = 3, ncol = 3), precision = prec, method = "slice")
  expect_equal(one, two)
  expect_true(is.matrix(one))
})

test_that("binomial slice sampling works for n=0", {
  m <- p <- rep_len(0, 10)
  k <- c(rep_len(0, 5), 1:5)
  n <- c(rep_len(0, 5), 2*(1:5))
  tau <- 10
  set.seed(20211203)
  norm <- rnorm(5, 0, sd = 1/sqrt(tau))
  set.seed(20211203)
  ss <- sample_binom_reg(p, k, n, m, tau, method = "slice")
  expect_equal(ss[1:5], norm)

  m <- p <- matrix(1, 4, 3)
  k <- matrix(c(rep_len(0, 6), 1:6), byrow = TRUE, nrow = 4)
  n <- 2*k
  tau <- diag(10, 3)
  set.seed(20211203)
  norm <- spam::rmvnorm.prec(2, rep(1, 3), Q = tau)
  set.seed(20211203)
  ss <- sample_binom_reg(p, k, n, m, tau, method = "slice")
  expect_equal(ss[1:2, ], norm)
})

# multinom -------------------------------------------------------------------

test_that("multivariate slice sampling works for multinomial", {
  tmp <- cbind(0, matrix(c(0, 0, 0.7, 0.7), nrow = 2))
  p <- array(rep(tmp, each = 400), dim = c(400, 2, 3))
  mean <- round(p)
  set.seed(99)
  n <- rpois(400*2, 1200)
  k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
  dim(k) <- c(400, 2, 3)
  stopifnot(rowSums(k, dims = 2) == n)

  z <- matrix(1, nrow = 2, ncol = 3)
  Q <- array(0, c(2, 2, 3))
  for(j in 1:3) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
  expect_error(sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q), "mean")
  mean <- mean[, , -1]
  expect_error(sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q), "precision")
  Q <- Q[, , -1]
  m <- apply(sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q), 3, mean)
  expect_true(all(abs(m - c(0, 0, log(2))) < 0.01))

  z <- cbind(1, c(0, 1), c(1, 1))
  expect_error(sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q))

  k[, 1, 2] <- 0
  m <- apply(ss <- sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q), 2:3, mean)
  expect_true(all(m[, 1] == 0))
  expect_true(all(abs(m[4:6] - c(0, log(2), log(2))) < 0.01))
})



test_that("multivariate slice sampling works 3-D z", {
  tmp <- cbind(0, matrix(c(0, 0, 0.7, 0.7), nrow = 2))
  p <- array(rep(tmp, each = 400), dim = c(400, 2, 3))
  mean <- round(p)
  set.seed(99)
  n <- rpois(400*2, 1200)
  k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
  dim(k) <- c(400, 2, 3)
  stopifnot(rowSums(k, dims = 2) == n)

  z <- matrix(1, nrow = 2, ncol = 3)
  Q <- array(0, c(2, 2, 3))
  for(j in 1:3) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
  mean <- mean[, , -1]
  Q <- Q[, , -1]
  set.seed(9909)
  m1 <- sample_multinom_reg(p = p, z = array(rep(z, each = 400), dim = dim(k)), k = k, mean = mean, precision = Q)
  set.seed(9909)
  m2 <- sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q)
  expect_equal(m1, m2)
})



test_that("multivariate slice sampling works when ref == 'last'", {
  tmp <- cbind(matrix(-0.7, nrow = 2, ncol = 2), 0)
  p <- array(rep(tmp, each = 400), dim = c(400, 2, 3))
  mean <- round(p[, , -1])
  set.seed(99)
  n <- rpois(400*2, 1200)
  k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
  dim(k) <- c(400, 2, 3)
  stopifnot(rowSums(k, dims = 2) == n)

  z <- matrix(1, nrow = 2, ncol = 3)
  Q <- array(0, c(2, 2, 2))
  for(j in 1:2) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
  m <- apply(sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q, ref = "last"), 3, mean)
  expect_true(all(abs(m - c(-log(2), -log(2), 0)) < 0.01))
})

test_that("multivariate slice sampling works when ref == {number}", {
  tmp <- cbind(matrix(-0.7, nrow = 2, ncol = 2), 0)
  p <- array(rep(tmp, each = 400), dim = c(400, 2, 3))
  mean <- round(p[, , -1])
  set.seed(99)
  n <- rpois(400*2, 1200)
  k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
  dim(k) <- c(400, 2, 3)
  stopifnot(rowSums(k, dims = 2) == n)

  z <- matrix(1, nrow = 2, ncol = 3)
  Q <- array(0, c(2, 2, 2))
  for(j in 1:2) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
  set.seed(20220124)
  m_last <- sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q, ref = "last")
  set.seed(20220124)
  m_2 <- sample_multinom_reg(p = p[, , c(1, 3, 2)], z = z, k = k[, , c(1, 3, 2)], mean = mean, precision = Q, ref = 2)
  set.seed(20220124)
  m_first <- sample_multinom_reg(p = p[, , c(3, 1:2)], z = z, k = k[, , c(3, 1:2)], mean = mean, precision = Q, ref = "first")

  expect_equal(m_2[, , c(2, 1, 3)], m_first)
  expect_equal(m_last[, , c(3, 1:2)], m_first)


})


test_that("multivariate slice sampling works when dimensions are 1", {
  k <- array(c(100, 1000, 10000), dim = c(1, 1, 3))
  p <- array(log(k / k[1]), dim = c(1, 1, 3))
  mean <- array(round(p[, , -1], 2), dim = c(1, 1, 2))
  set.seed(99)
  z <- matrix(1, nrow = 1, ncol = 3)
  Q <- array(1000, c(1, 1, 2))

  m <- sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q, ref = "first")
  expect_true(all(abs(m - p) < 0.01))
})












# other -------------------------------------------------------------------

test_that("One-dimensional precision matrices give the same results", {
  set.seed(20210119)
  one <- sample_binom_reg(0, 1, 100, mean = -3, precision = 1, method = "slice")
  set.seed(20210119)
  two <- sample_binom_reg(0, 1, 100, mean = -3, precision = matrix(1), method = "slice")
  expect_equal(one, two)


  set.seed(20210119)
  one <- sample_pois_reg(0, 10, mean = 3, precision = 1, method = "slice")
  set.seed(20210119)
  two <- sample_pois_reg(0, 10, mean = 3, precision = matrix(1), method = "slice")
  expect_equal(one, two)

})

