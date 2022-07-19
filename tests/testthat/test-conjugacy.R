test_that("Conjugacy works for 1-D cases", {
  mu1 <- conj_norm_mu(y = 1, tau = 10, mu0 = 0, tau0 = 1, params.only = TRUE)
  one <- matrix(1, 1, 1)
  ten <- matrix(10, 1, 1)
  mu2 <- conj_mvnorm_mu(y = 1, Q = ten, mu0 = 0, Q0 = one, params.only = TRUE)
  mu3 <- conj_matnorm_mu(y = one, V = ten, U = NULL, mu0 = 0, Q0 = one, params.only = TRUE)
  mu4 <- conj_matlm_beta(y = one, X = one, V = ten, U = NULL, mu0 = 0, Q0 = one, params.only = TRUE)
  mu5 <- conj_lm_beta(y = one, X = one, tau = 10, mu0 = 0, Q0 = one, params.only = TRUE)
  mu6 <- conj_matnorm_mu(y = one, V = one, U = ten, mu0 = 0, Q0 = one, params.only = TRUE)
  mu7 <- conj_diagmatlm_beta(y = one, X = one, V = ten, U = NULL, mu0 = 0, Q0 = one, params.only = TRUE)

  expect_equal(mu1$mu, mu2$mu)
  expect_equal(mu1$mu, mu3$mu)
  expect_equal(mu1$mu, mu4$mu)
  expect_equal(mu1$mu, mu5$mu)
  expect_equal(mu1$mu, mu6$mu)
  expect_equal(mu1$mu, mu7$mu)
})

test_that("Conjugacy works for 2-D independent cases", {
  mu1 <- c(
    conj_norm_mu(y = 1, tau = 10, mu0 = 0, tau0 = 1, params.only = TRUE)$mu,
    conj_norm_mu(y = 3, tau = 1, mu0 = 2, tau0 = 2, params.only = TRUE)$mu
  )
  one <- matrix(c(1, 3), 1, 2)
  two <- diag(c(1, 2))
  ten <- diag(c(10, 1))
  zero <- c(0, 2)
  X <- matrix(1, 1, 1)
  mu2 <- conj_mvnorm_mu(y = one, Q = ten, mu0 = zero, Q0 = two, params.only = TRUE)
  mu3 <- conj_matnorm_mu(y = one, V = ten, U = NULL, mu0 = zero, Q0 = two, params.only = TRUE)
  mu4 <- conj_matlm_beta(y = one, X = X, V = ten, U = NULL, mu0 = zero, Q0 = two, params.only = TRUE)
  mu6 <- conj_matnorm_mu(y = t(one), V = X, U = ten, mu0 = zero, Q0 = two, params.only = TRUE)
  mu7 <- conj_diagmatlm_beta(y = one, X = cbind(X, X), V = ten, U = NULL, mu0 = zero, Q0 = two, params.only = TRUE)

  expect_equal(mu1, mu2$mu)
  expect_equal(mu1, mu3$mu)
  expect_equal(mu1, mu4$mu)
  expect_equal(mu1, mu6$mu)
  expect_equal(mu1, mu7$mu)
})

test_that("diag = TRUE works for conj_matlm_beta", {
  V <- diag(1:3)
  X <- cbind(1, rep(c(0, 1), each = 5))
  set.seed(99)
  y <- matrix(rnorm(3*10), nrow = 10)
  mu0 <- matrix(0, nrow = 2, ncol = 3)
  Q0 <- diag(3) %x% diag(c(0.01, 0.1))
  set.seed(99)
  out <- conj_matlm_beta(y = y, X = X, V = V, mu0 = mu0, Q0 = Q0, diag = FALSE, use.chol = TRUE)
  set.seed(99)
  out2 <- conj_matlm_beta(y = y, X = X, V = V, mu0 = mu0, Q0 = Q0, diag = TRUE, use.chol = TRUE)
  set.seed(99)
  out3 <- conj_matlm_beta(y = y, X = X, V = V, mu0 = mu0, Q0 = Q0, diag = FALSE, use.chol = FALSE, params.only = TRUE)
  set.seed(99)
  out4 <- conj_matlm_beta(y = y, X = X, V = V, mu0 = mu0, Q0 = Q0, diag = TRUE, use.chol = FALSE, params.only = TRUE)

  expect_equal(out, out2)
  expect_equal(out3$mu, do.call(c, out4$mu))

  U <- diag(1:10)
  set.seed(99)
  out <- conj_matlm_beta(y = y, X = X, V = V, U = U, mu0 = mu0, Q0 = Q0, diag = FALSE, use.chol = TRUE)
  set.seed(99)
  out2 <- conj_matlm_beta(y = y, X = X, V = V, U = U, mu0 = mu0, Q0 = Q0, diag = TRUE, use.chol = TRUE)
  expect_equal(out, out2)


  set.seed(99)
  y <- matrix(rnorm(3*2), nrow = 2, ncol = 3)
  set.seed(99)
  out <- conj_matnorm_mu(y = y, V = V, mu0 = mu0, Q0 = Q0, diag = FALSE, use.chol = TRUE)
  set.seed(99)
  out2 <- conj_matnorm_mu(y = y, V = V, mu0 = mu0, Q0 = Q0, diag = TRUE, use.chol = TRUE)
  expect_equal(out, out2)

  U <- diag(c(1, 10))
  set.seed(99)
  out <- conj_matnorm_mu(y = y, V = V, U = U, mu0 = mu0, Q0 = Q0, diag = FALSE, use.chol = TRUE)
  set.seed(99)
  out2 <- conj_matnorm_mu(y = y, V = V, U = U, mu0 = mu0, Q0 = Q0, diag = TRUE, use.chol = TRUE)
  expect_equal(out, out2)
})



test_that("mu0 = NULL works for conj_matlm_beta", {
  V <- diag(1:3)
  X <- cbind(1, rep(c(0, 1), each = 5))
  set.seed(99)
  y <- matrix(rnorm(3*10), nrow = 10)
  mu0 <- matrix(0, nrow = 2, ncol = 3)
  Q0 <- diag(3) %x% diag(c(0.01, 0.1))

  out <- conj_matlm_beta(y = y, X = X, V = V, mu0 = mu0, Q0 = Q0, params.only = TRUE)
  out2 <- conj_matlm_beta(y = y, X = X, V = V, mu0 = NULL, Q0 = Q0, params.only = TRUE)
  out3 <- conj_matlm_beta(y = y, X = X, V = V, mu0 = mu0, Q0 = Q0, diag = TRUE, params.only = TRUE)
  out4 <- conj_matlm_beta(y = y, X = X, V = V, mu0 = NULL, Q0 = Q0, diag = TRUE, params.only = TRUE)

  expect_equal(out$mu, out2$mu)
  expect_equal(out$mu, do.call(c, out3$mu))
  expect_equal(out$mu, do.call(c, out4$mu))

  U <- diag(1:10)
  out <- conj_matlm_beta(y = y, X = X, V = V, U = U, mu0 = mu0, Q0 = Q0, params.only = TRUE)
  out2 <- conj_matlm_beta(y = y, X = X, V = V, U = U, mu0 = NULL, Q0 = Q0, params.only = TRUE)
  expect_equal(out, out2)


  set.seed(99)
  y <- matrix(rnorm(3*2), nrow = 2, ncol = 3)
  out <- conj_matnorm_mu(y = y, V = V, mu0 = mu0, Q0 = Q0, params.only = TRUE)
  out2 <- conj_matnorm_mu(y = y, V = V, mu0 = NULL, Q0 = Q0, params.only = TRUE)
  expect_equal(out, out2)

  out <- conj_mvnorm_mu(y, V, mu0 = rep_len(0, 3), Q0 = Q0[1:3, 1:3], params.only = TRUE)
  out2 <- conj_mvnorm_mu(y, V, mu0 = NULL, Q0 = Q0[1:3, 1:3], params.only = TRUE)
  expect_equal(out, out2)

})



test_that("newQ.chol= works", {
  Q <- cbind(1:2, c(2, 10))
  Q0 <- diag(1, 2)
  y <- cbind(1:3, 1)
  out1 <- conj_mvnorm_mu(y = y, Q = Q, Q0 = Q0, params.only = TRUE)
  out2 <- conj_mvnorm_mu(y = y, Q = Q, Q0 = Q0, newQ.chol = chol(3*Q + Q0), params.only = TRUE)
  expect_identical(out1, out2)
  set.seed(20220719)
  out1 <- conj_mvnorm_mu(y = y, Q = Q, Q0 = Q0)
  set.seed(20220719)
  out2 <- conj_mvnorm_mu(y = y, Q = Q, Q0 = Q0, newQ.chol = chol(3*Q + Q0))
  expect_equal(out1, out2)

  set.seed(20220719)
  out1 <- conj_mvnorm_mu(y = y, Q = Q, Q0 = Q0, mult = 1.25)
  set.seed(20220719)
  out2 <- conj_mvnorm_mu(y = y, Q = Q, Q0 = Q0, newQ.chol = chol(3*Q + 1.25*Q0))
  expect_equal(out1, out2)
})




