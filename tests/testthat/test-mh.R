prec <- matrix(c(10, 1, 2,
                 1, 20, 3,
                 2, 3, 30), byrow = TRUE, ncol = 3)

# binom -------------------------------------------------------------------


test_that("univariate MH sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- 10
  expect_error(mh_binom_reg(p, k, n, mean = 2, precision = 1), "'p' and 'k' and 'n' must all have")
  n <- rep(10, 3)
  set.seed(20210729)
  one <- mh_binom_reg(p, k, n, mean = 2, precision = 1, proposal_sd = 1)
  set.seed(20210729)
  two <- mh_binom_reg(p, k, n, mean = rep_len(2, 3), precision = 1, proposal_sd = 1)
  set.seed(20210729)
  three <- mh_binom_reg(p, k, n, mean = 2, precision = rep_len(1, 3), proposal_sd = 1)
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate MH sampling works for binomial", {
  p <- 1:3
  k <- c(3, 9, 10)
  n <- rep(10, 3)
  set.seed(20210729)
  one <- mh_binom_reg(p, k, n, mean = 2, precision = prec, proposal_sd = 1)
  set.seed(20210729)
  two <- mh_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec, proposal_sd = 1)
  set.seed(20210729)
  three <- mh_binom_reg(p, k, n, mean = 2, precision = prec, proposal_sd = 1)
  expect_equal(one, two)
  expect_equal(two, three)
  expect_true(!is.matrix(one))
})

test_that("multivariate MH sampling works for binomial (matrix)", {
  k <- c(3, 9, 10)
  p <- matrix(rep(1:3, each = 3), nrow = 3)
  n <- matrix(10, nrow = 3, ncol = 3)
  expect_error(mh_binom_reg(p, k, n, mean = 2, precision = prec, proposal_sd = 1),
               "'p' and 'k' and 'n' must all have")
  k <- matrix(rep(c(3, 9, 10), each = 3), nrow = 3)
  set.seed(20210729)
  one <- mh_binom_reg(p, k, n, mean = 2, precision = prec, proposal_sd = 1)
  set.seed(20210729)
  expect_error(mh_binom_reg(p, k, n, mean = rep_len(2, 3), precision = prec, proposal_sd = 1)
               , "'x' must be of length")
  two <- mh_binom_reg(p, k, n, mean = matrix(2, nrow = 3, ncol = 3), precision = prec, proposal_sd = 1)
  expect_equal(one, two)
  expect_true(is.matrix(one))
})



# other -------------------------------------------------------------------

test_that("One-dimensional precision matrices give the same results", {
  set.seed(20210119)
  one <- mh_binom_reg(0, 1, 100, mean = -3, precision = 1, proposal = "normal")
  set.seed(20210119)
  two <- mh_binom_reg(0, 1, 100, mean = -3, precision = matrix(1), proposal = "normal")
  expect_equal(one, two)


  set.seed(20210119)
  one <- ss_pois_reg(0, 10, mean = 3, precision = 1, proposal = "quad")
  set.seed(20210119)
  two <- ss_pois_reg(0, 10, mean = 3, precision = matrix(1), proposal = "quad")
  expect_equal(one, two)

})


test_that("n == 0 gives the same results for all methods", {
  p <- k <- n <- matrix(0, nrow = 5, ncol = 3)
  prec <- diag(3)
  set.seed(20210119)
  one <- mh_binom_reg(p, k, n, mean = -3, precision = prec, proposal = "normal")
  set.seed(20210119)
  two <- mh_binom_reg(p, k, n, mean = -3, precision = prec, proposal = "unif")
  set.seed(20210119)
  three <- mh_binom_reg(p, k, n, mean = -3, precision = prec, proposal = "quad")
  set.seed(20210119)
  four <- ss_binom_reg(p, k, n, mean = -3, precision = prec)
  expect_equal(one, two)
  expect_equal(one, three)
  expect_true(all(attr(one, "accept")))
  attr(one, "accept") <- NULL
  expect_equal(one, four)
})

test_that("Some n == 0 gives the same results for some methods", {

  p <- k <- n <- matrix(0:1, nrow = 1, ncol = 2)
  prec <- diag(100, 2)
  set.seed(20210120)
  three <- mh_binom_reg(p, k, n, mean = 2, precision = prec, proposal = "quad")
  set.seed(20210120)
  four <- ss_binom_reg(p, k, n, mean = 2, precision = prec)

  expect_true(attr(three, "accept")[1, 1])
  expect_equal(three[1, 1], four[1, 1])


  p <- k <- n <- matrix(1:0, nrow = 1, ncol = 2)
  prec <- diag(100, 2)
  set.seed(20210120)
  rnorm(1)
  three <- mh_binom_reg(p, k, n, mean = 2, precision = prec, proposal = "quad")
  set.seed(20210120)
  four <- ss_binom_reg(p, k, n, mean = 2, precision = prec)

  expect_true(attr(three, "accept")[1, 2])
  expect_equal(three[1, 2], four[1, 2])
})


