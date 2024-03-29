prec <- matrix(c(10, 1, 2,
                 1, 20, 3,
                 2, 3, 30), byrow = TRUE, ncol = 3)

# binom -------------------------------------------------------------------


test_that("univariate MH sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- 10
  expect_error(sample_binom_reg(p, k, n, mean = 2, precision = 1, method = "normal"), "'p' and 'k' and 'n' must all have")
  n <- rep(10, 3)
  set.seed(20210729)
  one <- sample_binom_reg(p, k, n, mean = 2, precision = 1, width = 1, method = "normal")
  set.seed(20210729)
  two <- sample_binom_reg(p, k, n, mean = rep_len(2, 3), precision = 1, width = 1, method = "normal")
  set.seed(20210729)
  three <- sample_binom_reg(p, k, n, mean = 2, precision = rep_len(1, 3), width = 1, method = "normal")
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate MH sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- rep(10, 3)
  set.seed(20210729)
  one <- sample_binom_reg(p, k, n, mean = 2, precision = prec, width = 1, method = "normal")
  set.seed(20210729)
  two <- sample_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec, width = 1, method = "normal")
  set.seed(20210729)
  three <- sample_binom_reg(p, k, n, mean = 2, precision = prec, width = 1, method = "normal")
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate MH sampling works for binomial (matrix)", {
  k <- c(3, 9, 10)
  p <- matrix(rep(1:3, each = 3), nrow = 3)
  n <- matrix(10, nrow = 3, ncol = 3)
  expect_error(sample_binom_reg(p, k, n, mean = 2, precision = prec, width = 1, method = "normal"),
               "'p' and 'k' and 'n' must all have")
  k <- matrix(rep(c(3, 9, 10), each = 3), nrow = 3)
  set.seed(20210729)
  one <- sample_binom_reg(p, k, n, mean = 2, precision = prec, width = 1, method = "normal")
  set.seed(20210729)
  expect_error(sample_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec, width = 1, method = "normal")
               , "'x' must be of length")
  two <- sample_binom_reg(p, k, n, mean = matrix(2, nrow = 3, ncol = 3), precision = prec, width = 1, method = "normal")
  expect_equal(one, two)
  expect_true(is.matrix(one))
})



# other -------------------------------------------------------------------

test_that("One-dimensional precision matrices give the same results", {
  set.seed(20210119)
  one <- sample_binom_reg(0, 1, 100, mean = -3, precision = 1, method = "normal")
  set.seed(20210119)
  two <- sample_binom_reg(0, 1, 100, mean = -3, precision = matrix(1), method = "normal")
  expect_equal(one, two)

  set.seed(20210119)
  one <- sample_binom_reg(0, 1, 100, mean = -3, precision = 1, method = "quad")
  set.seed(20210119)
  two <- sample_binom_reg(0, 1, 100, mean = -3, precision = matrix(1), method = "quad")
  set.seed(20210119)
  three <- sample_binom_reg(0, 1, 100, mean = -3, precision = matrix(1), method = "mv quad")
  expect_equal(one, two)
  expect_equal(one, three)


  set.seed(20210119)
  one <- sample_pois_reg(0, 1,  mean = -3, precision = 1, method = "normal")
  set.seed(20210119)
  two <- sample_pois_reg(0, 1, mean = -3, precision = matrix(1), method = "normal")
  expect_equal(one, two)

  set.seed(20210119)
  one <- sample_pois_reg(0, 1, mean = -3, precision = 1, method = "quad")
  set.seed(20210119)
  two <- sample_pois_reg(0, 1, mean = -3, precision = matrix(1), method = "quad")
  set.seed(20210119)
  three <- sample_pois_reg(0, 1, mean = -3, precision = matrix(1), method = "mv quad")
  expect_equal(one, two)
  expect_equal(one, three)

  set.seed(20210119)
  one <- sample_pois_reg(0, 1, mean = -3, precision = 1, method = "gamma")
  set.seed(20210119)
  two <- sample_pois_reg(0, 1, mean = -3, precision = matrix(1), method = "gamma")
  set.seed(20210119)
  three <- sample_pois_reg(0, 1, mean = -3, precision = matrix(1), method = "mv gamma")
  expect_equal(one, two)
  expect_equal(one, three)
})


test_that("n == 0 gives the same results for all methods", {
  p <- k <- n <- matrix(0, nrow = 5, ncol = 3)
  prec <- diag(3)
  set.seed(20210119)
  one <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "normal")
  set.seed(20210119)
  two <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "unif")
  set.seed(20210119)
  three <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "quad")
  set.seed(20210119)
  four <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "slice")
  set.seed(20210119)
  five <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "mv quad")
  set.seed(20210119)
  six <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "mv beta")
  set.seed(20210119)
  seven <- sample_binom_reg(p, k, n, mean = -3, precision = prec, method = "mv ind quad")
  expect_equal(one, two)
  expect_equal(one, three)
  expect_equal(one, five)
  expect_equal(one, six)
  expect_equal(one, seven)


  expect_true(all(attr(one, "accept")))
  attr(one, "accept") <- NULL

  expect_equal(one, four)
})

test_that("Some n == 0 gives the same results for some methods", {

  p <- k <- n <- matrix(0:1, nrow = 1, ncol = 2)
  prec <- diag(100, 2)
  set.seed(20210120)
  three <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "quad")
  set.seed(20210120)
  four <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "slice")

  expect_true(attr(three, "accept")[1, 1])
  expect_equal(three[1, 1], four[1, 1])


  p <- k <- n <- matrix(1:0, nrow = 1, ncol = 2)
  prec <- diag(100, 2)
  set.seed(20210120)
  rnorm(1)
  three <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "quad")
  set.seed(20210120)
  four <- sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "slice")

  expect_true(attr(three, "accept")[1, 2])
  expect_equal(three[1, 2], four[1, 2])
})


test_that("k == NA gives the same results for all methods", {
  k <- matrix(NA, nrow = 5, ncol = 3)
  L <- matrix(0, nrow = 5, ncol = 3)
  prec <- diag(3)
  set.seed(20210119)
  one <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "normal")
  set.seed(20210119)
  two <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "unif")
  set.seed(20210119)
  three <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "quad")
  set.seed(20210119)
  four <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "slice")
  set.seed(20210119)
  five <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "mv quad")
  set.seed(20210119)
  six <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "gamma")
  set.seed(20210119)
  seven <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "mv gamma")
  set.seed(20210119)
  eight <- sample_pois_reg(L, k, mean = -3, precision = prec, method = "slice", truncate = list(at = 4, allow = "below"))
  expect_equal(one, two)
  expect_equal(one, three)
  expect_equal(one, five)
  expect_equal(one, six)
  expect_equal(one, seven)

  expect_true(all(attr(one, "accept")))
  attr(one, "accept") <- NULL

  expect_equal(one, four)
  expect_equal(four, eight)

})

test_that("Some k == NA gives the same results for some methods", {

  k <- matrix(c(NA, 1), nrow = 1, ncol = 2)
  L <- matrix(0:1, nrow = 1, ncol = 2)
  prec <- diag(100, 2)
  set.seed(20210120)
  three <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "quad")
  set.seed(20210120)
  four <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "slice")
  set.seed(20210120)
  six <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "gamma")
  set.seed(20210120)
  eight <- sample_pois_reg(L, k, mean = 2, precision = prec, method = "slice", truncate = list(at = c(-1, 10), allow = "below"))

  expect_true(attr(three, "accept")[1, 1])
  expect_true(attr(six, "accept")[1, 1])
  expect_equal(three[1, 1], four[1, 1])
  expect_equal(three[1, 1], six[1, 1])
})


test_that("mean = -Inf works", {
  k <- c(0, NA, 0, 1)
  prec <- L <- c(1, 1, 1, 1)
  expect_error(sample_pois_reg(L, k, mean = -Inf, precision = prec, method = "normal"), "then k should be 0")
  expect_error(sample_pois_reg(L, k, mean = Inf, precision = prec, method = "normal"), "positive infinite")



  mn <- c(-Inf, -Inf, 0, 0)
  set.seed(20240118)
  one <- sample_pois_reg(L, k, mean = mn, precision = prec, method = "normal")
  set.seed(20240118)
  two <- sample_pois_reg(L, k, mean = mn, precision = prec, method = "unif")
  set.seed(20240118)
  three <- sample_pois_reg(L, k, mean = mn, precision = prec, method = "quad")
  set.seed(20240118)
  four <- sample_pois_reg(L, k, mean = mn, precision = prec, method = "slice")
  set.seed(20240118)
  six <- sample_pois_reg(L, k, mean = mn, precision = prec, method = "gamma")

  expect_error(sample_pois_reg(L, k, mean = mn, precision = matrix(prec, 2, 2), method = "mv quad"), "multivariate priors")
  expect_error(sample_pois_reg(L, k, mean = mn, precision = matrix(prec, 2, 2), method = "mv gamma"), "multivariate priors")
  expect_error(sample_pois_reg(L, k, mean = mn, precision = matrix(prec, 2, 2), method = "mv trunc"), "multivariate priors")
  expect_error(sample_pois_reg(L, k, mean = mn, precision = matrix(prec, 2, 2), method = "mv ind qu"), "multivariate priors")



  expect_equal(one[1:2], c(-Inf, -Inf))
  expect_equal(one[1:2], two[1:2])
  expect_equal(one[1:2], three[1:2])
  expect_equal(one[1:2], four[1:2])
  expect_equal(one[1:2], six[1:2])

})
