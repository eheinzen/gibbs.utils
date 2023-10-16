test_that("imputation works", {
  set.seed(99)
  try1 <- impute_conj_mvnorm_mu(
    y = c(1, 2),
    mu = c(1.1, 2.1),
    Q = matrix(c(1, 2, 2, 4), 2),
    mu0 = 1,
    tau0 = 0.1,
    impute = c(TRUE, FALSE)
  )
  set.seed(99)
  try2.param <- cond_mvnorm(
    y = c(1.1, 2.1),
    mu = c(1, 2),
    Q = matrix(c(1, 2, 2, 4), 2),
    which = 1,
    params.only = TRUE
  )
  try2 <- conj_norm_mu(
    y = try2.param$mu,
    tau = try2.param$tau,
    mu0 = 1,
    tau0 = 0.1
  )
  expect_equal(try1[1], try2)
})
